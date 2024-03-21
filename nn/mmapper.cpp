#include <iostream>
#include "../glib/misc/gregmmapper.hpp"

int main() {
    gtd::mmapper mapper{32'769};
    std::cout << mapper.capacity() << ", " << mapper.num_pages() << std::endl;
    size_t to_alloc = 16'385;
    mapper.reset(to_alloc);
    std::cout << "Allocated " << to_alloc << " bytes.\n";
    std::cout << "Num. pages = " << mapper.num_pages() << '\n';
    std::cout << "Capacity: " << mapper.capacity() << '\n';
    if (mapper)
        std::cout << "have memory" << std::endl;
    char *ptr = static_cast<char*>(mapper.get());
    uint64_t counter = to_alloc;
    srand(time(nullptr));
    while (counter --> 0)
        *ptr++ = 48 + rand() % 80;
    if (to_alloc)
        *--ptr = 0;
    printf("String: %s\n", mapper.get());
    mapper.reset();
    if (!mapper)
        printf("Mapper empty.\n");
    return 0;
}
