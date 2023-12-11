from supporting_functions import *

def create_whole_3dmol_text(chains):
    # Shuffled Brewer spectral scale from http://mkweb.bcgsc.ca/brewer/swatches/brewer.txt
    hexcolors = ['#9E0142', '#F46D43', '#FEE08B', '#ABDDA4', '#2B83BA', '#D53E4F', '#FC8D59', '#FFFFBF', '#99D594', '#3288BD', '#D7191C', '#FDAE61', '#E6F598', '#66C2A5', '#5E4FA2']
    text = "".join(("{\n",
        "\t\"rotate\":[{\"axis\":\"x\", \"angle\":-90},{\"axis\":\"y\", \"angle\":0},{\"axis\":\"z\", \"angle\":0}],\n",
        "\t\"sphere_transparency\": 0.45,\n",
        "\t\"chains\": [\n"))
    for nc in range(len(chains)):
        text += "".join(("\t\t{{\"chain\":\"{0}\",\n".format(chains[nc][-1]),
            "\t\t\t\"color\":\"{0}\",\n".format(hexcolors[nc%15]),
            "\t\t\t\"resid\":\"\"\n",
            "\t\t}"))
        if nc < len(chains) - 1:
            text += ",\n"
        else:
            text += "\n"
    text += "\t]\n"
    text += "}\n"
    return text


def create_chain_3dmol_text(chains, chain):
    # Shuffled Brewer spectral scale from http://mkweb.bcgsc.ca/brewer/swatches/brewer.txt
    hexcolors = ['#9e0142', '#f46d43', '#fee08b', '#abdda4', '#2b83ba', '#d53e4f', '#fc8d59', '#ffffbf', '#99d594', '#3288bd', '#d7191c', '#fdae61', '#e6f598', '#66c2a5', '#5e4fa2']

    text = "".join(("{\n",
        "\t\"rotate\":[{\"axis\":\"x\", \"angle\":-90},{\"axis\":\"y\", \"angle\":0},{\"axis\":\"z\", \"angle\":0}],\n",
        "\t\"sphere_transparency\": 0.45,\n",
        "\t\"chains\": [\n"))
    for nc in range(len(chains)):
        if chains[nc] == chain:
            text += "".join(("\t\t{{\"chain\":\"{0}\",\n".format(chains[nc][-1]),
                "\t\t\t\"color\":\"{0}\",\n".format(hexcolors[nc%15]),
                "\t\t\t\"resid\":\"\"\n",
                "\t\t}"))
            if nc < len(chains) - 1:
                text += ",\n"
            else:
                text += "\n"
        else:
            text += "".join(("\t\t{{\"chain\":\"{0}\",\n".format(chains[nc][-1]),
                "\t\t\t\"color\":\"{0}\",\n".format('#c8c8c8'),
                "\t\t\t\"transparency\":0.7,\n"
                "\t\t\t\"resid\":\"\"\n",
                "\t\t}"))
            if nc < len(chains) - 1:
                text += ",\n"
            else:
                text += "\n"
    text += "\t]\n"
    text += "}\n"
    return text


if __name__ == "__main__":
    from initialize_repository import *

    options, locations = initialize_repository()
    str_data = read_checkpoint(locations['FSYSPATH']['cache'] + 'str_data_completegen.pkl', locations['SYSFILES']['data_structure_template'])

    for pdbi in sorted(str_data):
        dmol_filename = locations['FSYSPATH']['3dmol_whole'] + pdbi + '.json'
        dmol_file = open(dmol_filename, 'w')
        dmol_file.write(create_whole_3dmol_text(sorted(str_data[pdbi]['ENCOMPASS']['structure']['kchains'])))
        dmol_file.close()

        for ch in str_data[pdbi]['ENCOMPASS']['structure']['kchains']:
            pdbi_ch = pdbi + '_' + ch
            dmol_filename = locations['FSYSPATH']['3dmol_chains'] + pdbi_ch + '.json'
            dmol_file = open(dmol_filename, 'w')
            dmol_file.write(create_chain_3dmol_text(sorted(str_data[pdbi]['ENCOMPASS']['structure']['kchains']), ch))
            dmol_file.close()
