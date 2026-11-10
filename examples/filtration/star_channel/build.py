from pore2process_demo.runtime.openlb import ROOT, build_case

BUILD = ROOT / '.local/star_channel/build'
BINARY = BUILD / 'capillary3d'

def main():
    return build_case(ROOT / 'examples/filtration/star_channel/openlb', BUILD, ROOT / 'results/star_channel/build')


if __name__ == '__main__':
    main()
