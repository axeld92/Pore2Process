from hydraulic.build import ROOT, build_case

BUILD = ROOT / '.local/star_channel/build'
BINARY = BUILD / 'capillary3d'

def main():
    return build_case(ROOT / 'star_channel/openlb', BUILD, ROOT / 'results/star_channel/build')


if __name__ == '__main__':
    main()
