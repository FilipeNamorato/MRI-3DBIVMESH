#!/bin/bash
#set -x

dir="MRI-3DBivMesh"

if [ $(basename "$PWD") = "$dir" ]; then
    if command -v git &> /dev/null; then
        echo "Cloning hexa-mesh-from-VTK repository..."
        git clone https://github.com/rsachetto/hexa-mesh-from-VTK.git
        if [ $? -eq 0 ]; then
            echo "Repository cloned successfully."
            cd hexa-mesh-from-VTK
            cmake .
            if [ $? -eq 0 ]; then
                echo "CMake configuration successful."
                make
                if [ $? -eq 0 ]; then
                    echo "Project build successful."
                else
                    echo "Failed to build project."
                    exit 1
                fi
            else
                echo "Failed to configure project with CMake."
                exit 1
            fi
        else
            echo "Failed to clone repository."
            exit 1
        fi
    else
        echo "Git is not installed. Unable to clone the repository."
        exit 1
    fi
else
    echo "You are not in the desired directory."
    echo "Please go to the ${dir} directory and run the script again."    
    exit 1
fi