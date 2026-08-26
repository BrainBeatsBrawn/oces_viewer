/*
 * Read OCES file and output data to stdout in compound ray .eye file format
 */

#include <iostream>
#include <string>

import oces.reader;
import oces.hexeye;

namespace local
{
    void stripFileSuffix (std::string& unixPath)
    {
        std::string::size_type pos (unixPath.rfind('.'));
        if (pos != std::string::npos) {
            // We have a '.' character
            std::string tmp (unixPath.substr (0, pos));
            if (!tmp.empty()) { unixPath = tmp; }
        }
    }
}

int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " path/to/oces_file.gltf\n";
        return -1;
    }
    std::string filename (argv[1]);

    oces::reader oces_reader (filename);
    oces::hexeye<float> heye (oces_reader.eye, 1.0f);
    local::stripFileSuffix (filename);
    heye.save (filename);

    return 0;
}
