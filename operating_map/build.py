from hydraulic.build import ROOT, build_case
BUILD = ROOT / '.local/operating_map/build_trt27'
BINARY = BUILD / 'capillary3d'
def main():
    return build_case(ROOT/'operating_map/openlb', BUILD, ROOT/'results/operating_map/build_trt27')
if __name__ == '__main__':
    main()
