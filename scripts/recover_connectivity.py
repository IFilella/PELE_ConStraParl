import argparse
import traceback

__author__ = "Carles Perez Lopez"


def parse_arguments():
    """
        Parse user arguments
        Output: list with all the user arguments
    """

    parser = argparse.ArgumentParser(description="""Description: Program to recover CONECT section of trajectories
     from prepared PDB files. """)
    required_named = parser.add_argument_group('required named arguments')
    required_named.add_argument("pdb_connected", type=str, help="PDB file with the CONECT section completed.")
    required_named.add_argument("--pdb_target",nargs="+",
                                help="List of PDB files without CONECT section, which must be added from "
                                     "the 'pdb_connected'."
                                )
    parser.add_argument("-lc", "--lig_chain_connected", type=str, default="L",
                        help="Ligand chain of the connected pdb. If it is not found, the program will select HETATM"
                             " without assigned chain.")
    parser.add_argument("-lt", "--lig_chain_target", type=str, default="L",
                        help="Ligand chain of the connected pdb. If it is not found, the program will select HETATM"
                             " without assigned chain.")
    parser.add_argument("-ov", "--overwrite", default=False, action='store_true',
                        help=" If set the input PDB targets will be overwritten by the connected version.")
    parser.add_argument("-rn", "--renum", default=False, action='store_true',
                        help=" If set the input PDB HETATOM indexes will be reassigned.")
    args = parser.parse_args()

    return args.pdb_connected, args.pdb_target, args.lig_chain_connected, args.lig_chain_target, args.overwrite, \
           args.renum


class LigandPDB:
    def __init__(self, pdb_file, ligand_chain="L"):
        self.pdb_file = pdb_file
        self.ligand_chain = ligand_chain
        self.content = self.read_pdb()

    def read_pdb(self):
        with open(self.pdb_file) as pdbf:
            content = pdbf.read()
        return content

    def get_models_dict(self):
        if "MODEL" in self.content:
            model_list = self.content.split("MODEL")
            model_dict = {}
            for n, model in enumerate(model_list[1:]):
                model_dict[n] = model
        else:
            model_dict = {}
            model = self.content
            model_dict[0] = model 
        return model_dict

    def get_model_ligand_content(self, model):
        lines = self.get_models_dict()[model].split("\n")
        ligand_lines = []
        ligand_blank = []
        for line in lines:
            if line.startswith("HETATM"):
                if line[21] == " ":
                    ligand_blank.append(line)
                if line[21] == self.ligand_chain or line[21] == " ":
                    ligand_lines.append(line)
        if ligand_blank:
            print("{}\nWARNING: We found that the ligand chain is a blank space. "
                  "Previuos lines will be included in the ligand!".format("\n".join(ligand_blank)))

        return "\n".join(ligand_lines)

    def get_connects(self):
        lines = self.content.split("\n")
        connect_lines = []
        for line in lines:
            if line.startswith("CONECT"):
                connect_lines.append(line)
        return "\n".join(connect_lines)

    def get_model_ligand_name_index_dictionary(self, model):
        ligand_lines = self.get_model_ligand_content(model=model).split("\n")
        lig_dict = {}
        for line in ligand_lines:
            lig_dict[line[12:16].strip().upper()] = line[6:11].strip()
        return lig_dict

    def add_connectivity(self, new_connects):
        if not self.get_connects():
            self.content = self.check_ending() + new_connects  # Update the content of the PDB
        else:
            print("Your PDB {} already contains a CONECT section!".format(self.pdb_file))
    
    def check_ending(self):
        lines = self.content.split("\n")
        if lines[-1].strip() == "END":
            print("END detected")
            del lines[-1]
            return "\n".join(lines)
        elif lines[-2].strip() == "END":
            print("END detected")
            del lines[-2]
            return "\n".join(lines)
        else:
            return self.content


def renum_hetatom(pdb_content):
    lines = pdb_content.split("\n")
    new_pdb = []
    for line in lines:
        if line.startswith("ATOM  ") or line.startswith("TER  "):
            try:
                index = int(line[6:11])
            except IndexError:
                print("TER not found in the last atom.")
        if line.startswith("HETATM"):
            curr_index = int(line[6:11])
            new_index = index + curr_index
            new_line = list(line)
            new_line[6:11] = f"{new_index:>5}"
            line = "".join(new_line)
        new_pdb.append(line)
    return "\n".join(new_pdb)


def create_index_relation_between_connect_and_to_connect(pdb_connected, pdb2connect):
    pdb_connected_names_dict = pdb_connected.get_model_ligand_name_index_dictionary(0)
    pdb2_connect_names_dict = pdb2connect.get_model_ligand_name_index_dictionary(0)
    index_relations = {}
    for key, value in pdb2_connect_names_dict.items():
        index_relations[pdb_connected_names_dict[key]] = value
    return index_relations


def rebuild_connect_line(connections_list):
    new_connect_pattern = "{:5s}"*len(connections_list)
    new_connect_line = "CONECT " + new_connect_pattern.format(*connections_list)
    return new_connect_line


def recover_connectivity(pdb_connected, pdb_to_connect):
    index_relations = create_index_relation_between_connect_and_to_connect(pdb_connected, pdb_to_connect)
    connects_pdb_connected = pdb_connected.get_connects().split("\n")
    connects_lines = []
    for line in connects_pdb_connected:
        connections_list = line.split()[1:]
        for n, index in enumerate(connections_list):
            try:
                new_index = index_relations[index]
                connections_list[n] = new_index
            except KeyError:
                print("Index {} not found in Ligand".format(index))
        new_connect = rebuild_connect_line(connections_list)
        connects_lines.append(new_connect)
    connects = "\n".join(connects_lines)
    return connects


def main(pdb_connected, pdb_to_connect, ligand_chain_connected="L", ligand_chain_to_connect="L", overwrite=False,
         renum=False):
    pdbc = LigandPDB(pdb_file=pdb_connected, ligand_chain=ligand_chain_connected)
    print("Ligand: {}".format(pdb_connected))
    for p2c in pdb_to_connect:
        pdb2c = LigandPDB(pdb_file=p2c, ligand_chain=ligand_chain_to_connect)
        if renum:
            old_pdb = pdb2c.content
            pdb2c.content = renum_hetatom(old_pdb)
        new_connects = recover_connectivity(pdbc, pdb2c)
        pdb2c.add_connectivity(new_connects)
        if overwrite:
            with open(p2c, "w") as out:
                out.write(pdb2c.content)
        else:
            if "/" in p2c:
                path = p2c.split("/")
                path[-1] = "conn"+path[-1]
                new_p2c = "/".join(path)
            else:
                new_p2c = p2c
            with open(new_p2c, "w") as out:
                out.write(pdb2c.content)

if __name__ == '__main__':
    pdbc, pdbt, lig_ch_c, lig_ch_t, ov, renum = parse_arguments()
    main(pdbc, pdbt, lig_ch_c, lig_ch_t, ov, renum)


