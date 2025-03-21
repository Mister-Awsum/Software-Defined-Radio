/*
Comp Eng 3DY4 (Computer Systems Integration Project)

Department of Electrical and Computer Engineering
McMaster University
Ontario, Canada
*/

#include "iofunc.h"

// Some basic functions for printing information from vectors
// or to read from/write to binary files in real format

void printRealVector(const std::vector<dy4::real> &x) {
	std::cerr << "Printing real vector of size " << x.size() << "\n";
	for (int i = 0; i < (int)x.size(); i++) {
		std::cerr << x[i] << " ";
	}
	std::cerr << "\n";
}

void printComplexVector(const std::vector<std::complex<dy4::real>> &X) {
	std::cerr << "Printing complex vector of size " << X.size() << "\n";
	for (int i = 0; i < (int)X.size(); i++) {
		std::cerr << X[i] << "";
	}
	std::cerr << "\n";
}

// Assumes data in the raw binary file is in the format corresponding to 'real'
void readBinData(const std::string in_fname, std::vector<dy4::real> &bin_data)
{
	std::ifstream fdin(in_fname, std::ios::binary);
	if (!fdin) {
		std::cerr << "File " << in_fname << " not found ... exiting\n";
		exit(1);
	} else {
		std::cerr << "Reading raw binary from \"" << in_fname << "\"\n";
	}

	fdin.seekg(0, std::ios::end);
	const unsigned int num_samples = fdin.tellg() / sizeof(dy4::real);

	bin_data.resize(num_samples);
	fdin.seekg(0, std::ios::beg);
	fdin.read(reinterpret_cast<char *>(&bin_data[0]), num_samples * sizeof(dy4::real));
	fdin.close();
}

// Assumes data in the raw binary file is in the format corresponding to 'real'
void writeBinData(const std::string out_fname, const std::vector<dy4::real> &bin_data)
{
	// std::cerr << "Writing raw binary to \"" << out_fname << "\"\n";
	std::ofstream fdout(out_fname, std::ios::binary);
	for (int i = 0; i < (int)bin_data.size(); i++) {
		fdout.write(reinterpret_cast<const char *>(&bin_data[i]), sizeof(bin_data[i]));
	}
	fdout.close();
}

//////////////////////////////////////////////////////////////
// New code as part of benchmarking/testing and the project
//////////////////////////////////////////////////////////////

void generate_random_values(std::vector<dy4::real>& x, const dy4::real& lower_bound = -1., const dy4::real& upper_bound = 1.) {

	for (int i = 0; i < (int)x.size(); i++) {
		x[i] = lower_bound + (upper_bound - lower_bound) * ((dy4::real)(rand()) / RAND_MAX);
	}
}