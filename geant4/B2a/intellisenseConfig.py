import os
import json

results = []

def main():
    """
    This function finds the folders with header files and source files
    """
    head = '/home/qifan/projects/partiTrans/geant4-v11.1.1/source'
    pop(head)
    # print(results)

    cpp_config = {
    "configurations": [
        {
            "name": "Linux",
            "includePath": [
                "${workspaceFolder}/**"
            ],
            "defines": [],
            "compilerPath": "/usr/bin/gcc",
            "cStandard": "c17",
            "cppStandard": "gnu++14",
            "intelliSenseMode": "linux-gcc-x64",
            "configurationProvider": "ms-vscode.makefile-tools"
        }
    ],
    "version": 4
    }
    cpp_config['configurations'][0]['includePath'].extend(results)
    targetFile = '.vscode/c_cpp_properties.json'
    json_object = json.dumps(cpp_config, indent=4)
    with open(targetFile, 'w') as f:
        f.write(json_object)


def pop(folder):
    """
    This is the recursion funciton to populuate results
    """
    global results
    subs = os.listdir(folder)
    for sub in subs:
        if '.hh' in sub or '.cc' in sub:
            results.append(folder)
            break
    for sub in subs:
        sub_full = os.path.join(folder, sub)
        if os.path.isdir(sub_full):
            pop(sub_full)


if __name__ == '__main__':
    main()