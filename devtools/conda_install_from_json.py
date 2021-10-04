import sys
import json
import delegator
import os


def read_requirements(fnames):
    keys = ["channels", "packages", "cpu_packages", "gpu_packages"]
    with open(fnames[0]) as fin:
        d = json.load(fin)
    for fname in fnames[1:]:
        with open(fname) as fin:
            d2 = json.load(fin)
            for key in keys:
                d[key].extend(d2[key])
    return d


def add_platform_specific_packages(cmd, requirements):
    if 'CPU_ONLY' in os.environ and os.environ['CPU_ONLY'] == '1':
        for package in requirements['cpu_packages']:
            cmd = '%s "%s"' % (cmd, package)
        return cmd
    for package in requirements['gpu_packages']:
        cmd = '%s "%s"' % (cmd, package)
    return cmd


def add_wheels(requirements):
    wheels = requirements.get('wheels', [])
    if 'CPU_ONLY' in os.environ and os.environ['CPU_ONLY'] == '1':
        wheels.extend(requirements.get('cpu_wheels', []))
    else:
        wheels.extend(requirements.get('gpu_wheels', []))
    for whl in wheels:
        cmd = "pip install %s" % whl
        c = delegator.run(cmd)
        print(c.out)
        print(c.err)


def main(fnames):
    requirements = read_requirements(fnames)
    cmd = "conda install -y -q "
    for channel in requirements['channels']:
        cmd = "%s -c %s" % (cmd, channel)
    for package in requirements['packages']:
        cmd = '%s "%s"' % (cmd, package)
    cmd = add_platform_specific_packages(cmd, requirements)
    print("Conda Command")
    print("%s" % cmd)
    c = delegator.run(cmd)
    print(c.out)
    print(c.err)
    add_wheels(requirements)
    sys.exit(c.return_code)


if __name__ == "__main__":
    main(sys.argv[1:])
