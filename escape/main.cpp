#include "../glib/nbod/gregsys.hpp"
#include "../glib/nbod/gregastro.hpp"
#include <iomanip>

/* Test program to test escape velocity formulae that I've derived myself. */

// #define G 0.0000000000667430l

#define DELTA_T (0.125l/65'536.0l)
#define ITERATIONS 1'000'000

long double escape_speed_v1_com(long double mass1, long double mass2, long double r) {
	return mass2*sqrtl((2*gtd::sys::G_SI)/(r*(mass1 + mass2)));
}

long double escape_speed_v2_com(long double mass1, long double mass2, long double r) {
	return mass1*sqrtl((2*gtd::sys::G_SI)/(r*(mass1 + mass2)));
}

int main() {
	long double m1 = 1'000'000.0l;
	long double m2 = 2'000'000.0l;
	long double r = 4;
	long double v1_com = escape_speed_v1_com(m1, m2, r);
	long double v2_com = escape_speed_v2_com(m1, m2, r); /*
	gtd::bod_0f b1{m1, 0, {-(r*m2)/(m1 + m2), 0, 0}, {-v1_com, 0, 0}};
	gtd::bod_0f b2{m2, 0, {(r*m1)/(m1 + m2), 0, 0}, {v2_com, 0, 0}}; */
    std::mt19937_64 rng{std::random_device{}()};
    std::uniform_real_distribution<long double> dist{0, 1};
    gtd::vec3 v1{dist(rng), dist(rng), dist(rng)};
    gtd::vec3 v2 = -v1; // must be in opposite directions to keep COM vel at 0
    v1 *= v1_com/v1.magnitude();
    v2 *= v2_com/v2.magnitude();
    std::cout << "v1_com: " << v1_com << ", v1.magnitude(): " << v1.magnitude() << std::endl;
    std::cout << "v2_com: " << v2_com << ", v2.magnitude(): " << v2.magnitude() << std::endl;
    std::cout << "v1.unit_vector(): " << v1.unit_vector() << std::endl;
    std::cout << "v2.unit_vector(): " << v2.unit_vector() << std::endl;
    gtd::bod_0f b1{m1, cbrtl(m1/m2), {-(r*m2)/(m1 + m2), 0, 0}, v1};
    gtd::bod_0f b2{m2, 1, {(r*m1)/(m1 + m2), 0, 0}, v2};
	gtd::sys sys;
	sys.add_body(b1);
	sys.add_body(b2);
	sys.timestep(DELTA_T);
	sys.iters(ITERATIONS);
    std::cout << "sys.com_pos(): " << sys.com_pos() << std::endl;
    std::cout << "sys.com_vel(): " << sys.com_vel() << std::endl;
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
    drdt << "time(s),dr/dt(m/s)\r\n";
    long double two_h = 2*h;
    gtd::vec3 val;
    gtd::image_dimensions dims{2000, 2000};
    gtd::asc_0f scene{dims.x, dims.y};
    scene.add_body(sys[0]);
    scene.add_body(sys[1]);
    std::cout << sys[0] << "\n---------------\n" << sys[1] << std::endl;
    gtd::star_t s1{0, 0, {0, 0, 1'000}, {}, 1, 1};
    scene.add_star(s1);
    gtd::cam cam;
    cam.set_position({0, 0, 50});
    cam.set_direction({0, 0, -1});
    cam.set_image_dimensions(dims);
    scene.follow_camera(&cam);
    scene.set_num_decor_stars(0);
    std::ostringstream oss;
    unsigned int num_zeros = ((unsigned int) ceill(logl(epochs)));
    scene.render();
    oss.fill(48);
    oss << "epoch" << std::setw(num_zeros) << counter << ".bmp";
    std::string path = oss.str();
    scene.write(path.c_str());
	while (counter < epochs) {
		std::cout << "Epoch " << counter << "\n------------" << std::endl;
		sys.evolve();
		std::cout << "Positions:\n";
		for (const auto &bod : sys)
			std::cout << (val = bod.pos()) << '\n' << val.magnitude() << std::endl;
		std::cout << "Velocities:\n";
		for (const auto &bod : sys)
			std::cout << (val = bod.vel()) << '\n' << val.magnitude() << std::endl;
        std::cout << "Accelerations:\n";
        for (const auto &bod : sys)
            std::cout << (val = bod.acceleration()) << '\n' << val.magnitude() << std::endl;
		std::cout << "Separation: " << (sep = (sys[0].pos() - sys[1].pos()).magnitude()) << std::endl;
		separations.push_back(sep);
		if (!counter)
			sep_drdt.push_back((separations[1] - separations[0])/(h));
		else if (counter < (epochs - 1))
			sep_drdt.push_back((separations[counter + 1] - separations[counter - 1])/(two_h));
		else
			sep_drdt.push_back((separations.back() - separations[counter])/(h));
        drdt << counter++*h << ',' << sep_drdt.back() << "\r\n";
        if (!(counter % 10)) {
            oss.seekp(5, std::ios_base::beg);
            oss << std::setw(num_zeros) << counter;
            path = oss.str();
            scene.render();
            scene.write(path.c_str());
        }
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
