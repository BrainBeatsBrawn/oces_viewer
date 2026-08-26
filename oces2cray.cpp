/*
 * Read OCES file and output data to stdout in compound ray .eye file format
 */

#include <iostream>
#include <string>

import oces.reader;

int main (int argc, char** argv)
{
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " path/to/oces_file.gltf\n";
        return -1;
    }
    std::string filename (argv[1]);

    oces::reader oces_reader (filename);
    oces_reader.eye.output_compound_ray_csv();

    return 0;
}
