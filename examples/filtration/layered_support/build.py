"""Build the support reference outside the original OpenLB distribution."""
from pore2process_demo.runtime.openlb import build_case, ROOT

BUILD = ROOT / ".local/filtration/build"
BINARY = BUILD / "capillary3d"


def main():
    return build_case(ROOT / "examples/filtration/layered_support/openlb", BUILD, ROOT / "results/filtration")


if __name__ == "__main__":
    main()
