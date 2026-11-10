"""Build only project-owned files; verify the original OpenLB archive."""
from pore2process_demo.runtime.openlb import ROOT, build_case

BUILD=ROOT/".local/graded_support/build"
BINARY=BUILD/"capillary3d"

def main():
    return build_case(ROOT/"examples/filtration/graded_support/openlb",BUILD,ROOT/"results/graded_support")

if __name__=="__main__":
    main()
