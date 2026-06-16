from conan import ConanFile
from conan.tools.cmake import CMake, CMakeToolchain, CMakeDeps, cmake_layout


class PolyhedralTemplateMatchingConan(ConanFile):
    name = "polyhedral-template-matching"
    version = "2.0.2"
    package_type = "static-library"
    license = "MIT"
    settings = "os", "arch", "compiler", "build_type"
    default_options = {"hwloc/*:shared": True}
    requires = (
        "boost/1.88.0",
        "onetbb/2021.12.0",
        "common-neighbor-analysis/[>=2.0]",
        "coretoolkit/[>=2.5]",
        "structure-identification/[>=2.1]",
        "spdlog/1.14.1",
        "nlohmann_json/3.11.3",
        "yaml-cpp/0.8.0",
    )
    exports_sources = "CMakeLists.txt", "include/*", "src/*", "dependencies/*"

    def layout(self):
        cmake_layout(self)

    def generate(self):
        CMakeToolchain(self).generate()
        CMakeDeps(self).generate()

    def build(self):
        cmake = CMake(self)
        cmake.configure()
        cmake.build()

    def package(self):
        cmake = CMake(self)
        cmake.install()

    def package_info(self):
        self.cpp_info.set_property(
            "cmake_target_name",
            "polyhedral-template-matching::polyhedral-template-matching",
        )
        self.cpp_info.libs = ["polyhedral-template-matching_lib"]
        self.cpp_info.requires = [
            "boost::headers",
            "onetbb::onetbb",
            "common-neighbor-analysis::common-neighbor-analysis",
            "coretoolkit::coretoolkit",
            "structure-identification::structure-identification",
            "nlohmann_json::nlohmann_json",
            "spdlog::spdlog",
            "yaml-cpp::yaml-cpp",
        ]
