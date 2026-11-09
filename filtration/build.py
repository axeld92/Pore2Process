"""Build the support reference outside the original OpenLB distribution."""
from hydraulic.build import build_case, ROOT

BUILD = ROOT / ".local/filtration/build"
BINARY = BUILD / "capillary3d"


def main():
    return build_case(ROOT / "filtration/openlb", BUILD, ROOT / "results/filtration")


if __name__ == "__main__":
    main()
