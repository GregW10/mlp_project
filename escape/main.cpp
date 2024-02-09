#include "../glib/nbod/gregsys.hpp"

/* Throw-away program to test escape velocities */

#define G 0.0000000000667430l

#define DELTA_T 0.125l/65'536.0l
#define ITERATIONS 1'000'000

long double escape_speed_v1_com(long double mass1, long double mass2, long double r) {
	return mass2*sqrtl((2*G)/(r*(mass1 + mass2)));
}

long double escape_speed_v2_com(long double mass1, long double mass2, long double r) {
	return mass1*sqrtl((2*G)/(r*(mass1 + mass2)));
}

int main() {
	long double m1 = 1'000'000.0l;
	long double m2 = 2'000'000.0l;
	long double r = 2;
	long double v1_com = escape_speed_v1_com(m1, m2, r);
	long double v2_com = escape_speed_v2_com(m1, m2, r);
	gtd::bod_0f b1{m1, 0, {-(r*m2)/(m1 + m2), 0, 0}, {-v1_com, 0, 0}};
	gtd::bod_0f b2{m2, 0, {(r*m1)/(m1 + m2), 0, 0}, {v2_com, 0, 0}};
	gtd::sys sys;
	sys.add_body(b1);
	sys.add_body(b2);
	sys.set_timestep(DELTA_T);
	sys.set_iterations(ITERATIONS);
	unsigned int epochs = 100'000;
    long double h = DELTA_T*ITERATIONS; // epoch timestep
    std::cout << "Time-step: " << DELTA_T << " s\nIterations per epoch: " << ITERATIONS << "\nNumber of epochs: " <<
    epochs << "\nEpoch time-step: " << h << " s\n------------------------\n" << std::endl;
	unsigned int counter = 0;
	std::cout << "Initial Positions:\n";
	for (const auto &bod : sys)
		std::cout << bod.pos() << std::endl;
	std::cout << "Initial Velocities:\n";
	for (const auto &bod : sys)
		std::cout << bod.vel() << std::endl;
	long double sep;
	std::cout << "Initial Separation: " << (sep = (sys[0].pos() - sys[1].pos()).magnitude()) << '\n' << std::endl;
	std::vector<long double> separations;
	separations.reserve(epochs + 1);
	separations.push_back(sep);
	std::vector<long double> sep_drdt;
	sep_drdt.reserve(epochs - 1);
    std::ofstream drdt{"drdt_approx.csv", std::ios_base::out | std::ios_base::trunc};
    drdt << "time,dr/dt\r\n";
    long double two_h = 2*h;
	while (counter < epochs) {
		std::cout << "Epoch " << counter << "\n------------" << std::endl;
		sys.evolve();
		std::cout << "Positions:\n";
		for (const auto &bod : sys)
			std::cout << bod.pos() << std::endl;
		std::cout << "Velocities:\n";
		for (const auto &bod : sys)
			std::cout << bod.vel() << std::endl;
        std::cout << "Accelerations:\n";
        for (const auto &bod : sys)
            std::cout << bod.acceleration() << std::endl;
		std::cout << "Separation: " << (sep = (sys[0].pos() - sys[1].pos()).magnitude()) << std::endl;
		separations.push_back(sep);
		if (!counter)
			sep_drdt.push_back((separations[1] - separations[0])/(h));
		else if (counter < (epochs - 1))
			sep_drdt.push_back((separations[counter + 1] - separations[counter - 1])/(two_h));
		else
			sep_drdt.push_back((separations.back() - separations[counter])/(h));
        drdt << counter++*h << ',' << sep_drdt.back() << "\r\n";
		putchar('\n');
	}
    drdt.close();
    drdt.open("drdt_approx.bin", std::ios_base::out | std::ios_base::trunc | std::ios_base::binary);
    uint64_t size = sizeof(long double);
    drdt.write((char *) &size, sizeof(uint64_t));
    drdt.write((char *) &h, sizeof(long double));
    drdt.write((char *) sep_drdt.data(), sizeof(long double)*sep_drdt.size());
    drdt.close();
	return 0;
}
