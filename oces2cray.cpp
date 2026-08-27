/*
 * Read OCES file and output data to stdout in compound ray .eye file format
 */

#include <iostream>
#include <string>

import oces.reader;

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
        std::cout << "Usage: " << argv[0] << " path/to/oces_file.gltf [-6]\n--\n"
                  << "Reads OCES file and writes the compound-ray .eye file.\n"
                  << "Also optionally writes a compound-ray file for the hex-equivalent\n"
                  << "eye, along with an sm::hexgrid HDF5 file if the -6 option is passed.\n";
        return -1;
    }
    std::string filename (argv[1]);

    oces::reader oces_reader (filename);

    std::string filebase = filename;
    local::stripFileSuffix (filebase);

    std::string fileeye = filebase + ".eye";
    oces_reader.eye.output_compound_ray_csv (fileeye);

    if (argc > 2) {
        std::string harg = std::string(argv[2]);
        if (harg == std::string("-6")) {
            oces_reader.setup_hexeye (1.0f);
            oces_reader.heye.save (filebase);
        }
    }

    return 0;
}
