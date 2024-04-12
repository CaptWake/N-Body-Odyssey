#include <iostream>

int main(int argc, char **argv) {
    if (argc != 2) {
        std::cout << "Usage: " << argv[0] << " <number>" << std::endl;
        return 1;
    }
    std::cout << "new version release" << std::endl;
    for (int i = 0; i < std::atoi(argv[1]); i++) {
        //execute scpd-project new version
        std::string command = "\"/home/dadi/Documents/Scuola/Informatica/Magistrale/Anno 1/Semestre 2/Sistemi paralleli e distribuiti/Progetto/SCPD-Project/cmake-build/release/SCPD_Project\"";
        std::string param = "-c \"/home/dadi/Documents/Scuola/Informatica/Magistrale/Anno 1/Semestre 2/Sistemi paralleli e distribuiti/Progetto/asd.json\"";
        command += " " + param;
        std::system(command.c_str());
    }
    return 0;
}