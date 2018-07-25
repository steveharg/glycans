import sys, json
import matplotlib.pyplot as plt

# These enumeration-type variables are used to gather metrics on the reasons why some wurcs strings are rejected
WURCS2_STR_IDX = 0
UNCERTAIN_LIN_COUNT_IDX = 1
STATISTIC_LIP_IDX = 2
MORE_THAN_TWO_SIDED_LINK_IDX = 3
ALTERNATIVE_LIP_IDX = 4
UNCERTAIN_CARBON_NUM_IDX = 5
SINGLE_SIDED_LINK_IDX = 6
OTHER_LINK_PROBLEM_IDX = 7

def parse_glytoucan_results(json_file, num_datafile_lines):

    # Get glycan sequences from json file (previously generated via sparql query in glytoucan_rdf_to_json.py)
    with open(json_file) as f:
        results = json.load(f)

    if num_datafile_lines != None:
        # Reduce the dataset size whilst testing
        # results["results"]["bindings"] = results["results"]["bindings"][6000:6020]
        results["results"]["bindings"] = results["results"]["bindings"][0:num_datafile_lines-1]

    single_sided_link_count = 0
    double_sided_link_count = 0
    total_link_count = 0

    # the total number of reactions (child+link+parent) found across all glycans
    reaction_count = 0

    # dictionary of structures (key = PrimaryId) and their reactions (value = list of reactions present)
    reactions_bag = {}

    # dictionary of structures (key = PrimaryId) which we can't use, and flags indicating reasons why
    # [uncertain_lin_count]
    rejected_structures = {}

    # parse the results
    for result in results["results"]["bindings"]:
        statistic_lip = False
        alternative_lip = False
        uncertain_lin_count = False

        # The full wurcs v2 string, in this format: <Version>/<Unit Count>/<UniqueRES List>/<RES Sequence>/<LIN List>
        wurcs_string = result["Sequence"]["value"]

        # check that the '£' character isn't used anywhere, because we'll use it as a separator later
        if '£' in wurcs_string:
            print("£ in wurcs string")
            sys.exit(0)

        primary_id = result["PrimaryId"]["value"]

        wurcs_string_split= wurcs_string.split('/', 2)
        wurcs2_header = wurcs_string_split[0]
        wurcs2_unit_count = wurcs_string_split[1]
        wurcs2_unit_count_split = wurcs2_unit_count.split(',')
        num_unique_residues = int(wurcs2_unit_count_split[0])
        num_residues = int(wurcs2_unit_count_split[1])

        num_links_str = wurcs2_unit_count_split[2]

        # If number of links is uncertain, reject
        if num_links_str[-1] == '+':
            num_links_str = num_links_str[:-1]
            uncertain_lin_count = True
            if primary_id not in rejected_structures:
                rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

            rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
            rejected_structures[primary_id][UNCERTAIN_LIN_COUNT_IDX] = uncertain_lin_count

        num_links = int(num_links_str)

        # The whole wurcs string, apart from <Version>/<Unit Count>/
        # i.e. just: <UniqueRES List>/<RES Sequence>/<LIN List>
        wurcs2_res_list_res_seq_lin_list = wurcs_string_split[2]
        wurcs2_res_list_res_seq_lin_list_split = wurcs2_res_list_res_seq_lin_list.split(']')

        # Make a list of the unique residues
        wurcs2_unique_res_list = []
        for unique_res in wurcs2_res_list_res_seq_lin_list_split:
            wurcs2_unique_res_list.append(unique_res[1:])

        # The last item in this list will be <RES Sequence>/<LIN List> which we don't want here, but do want elsewhere,
        # so pop into a new variable
        wurcs2_res_seq_lin_list = wurcs2_unique_res_list.pop()

        # Split <RES Sequence>/<LIN List> by /
        wurcs2_res_seq_lin_list_split = wurcs2_res_seq_lin_list.split('/', 1)
        # Get the <RES Sequence>
        wurcs2_res_seq_str = wurcs2_res_seq_lin_list_split[0]

        # <RES Sequence> is formed of hyphon-separated numbers, i.e. <RES #1>-<RES #2>-...-<RES #n>
        # So split by - to get the list of residue sequence numbers
        wurcs2_res_seq_str_split = wurcs2_res_seq_str.split('-')
        wurcs2_res_seq_list = []
        for res_seq_str in wurcs2_res_seq_str_split:
            wurcs2_res_seq_list.append(int(res_seq_str))

        # Use the residue sequence number list to build a list of actual residues
        wurcs2_res_list = []
        for res in wurcs2_res_seq_list:
            wurcs2_res_list.append(wurcs2_unique_res_list[res-1])

        if num_links > 0:

            wurcs2_lin_list = wurcs2_res_seq_lin_list_split[1]

            # wurcs2_res_seq_lin_list_split should contain two items (res_seq and lin_list) - if there are more, exit
            if len(wurcs2_res_seq_lin_list_split) > 2:
                print("more than two items in wurcs2_res_seq_lin_list_split: ", wurcs2_res_seq_lin_list_split)
                print("exiting")
                sys.exit(0)

            # Each link is separated by an undrscore
            wurcs2_lin_list_split = wurcs2_lin_list.split('_')
            # Create a list of the links
            link_list = []
            for link in wurcs2_lin_list_split:

                # If links only have a percentage certainty, reject
                if '%' in link:
                    statistic_lip = True
                    if primary_id not in rejected_structures:
                        rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                    rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                    rejected_structures[primary_id][STATISTIC_LIP_IDX] = statistic_lip

                total_link_count += 1
                link_list.append(link)
                # The two ends (child and parent) of the link are separated by a dash
                link_split = link.split('-')
                # There must be at least two ends to the link
                if len(link_split) > 1:
                    # If there are more than two 'ends' (i.e. multiple parents or children)
                    # that's valid, but we ignore for now
                    if len(link_split) > 2:
                        # print("warning - more than two elements in link")

                        if primary_id not in rejected_structures:
                            rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                        rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                        rejected_structures[primary_id][MORE_THAN_TWO_SIDED_LINK_IDX] = True
                    # Reject cases where the are multiple alternative links (left side)
                    if '|' in link_split[0]:
                        # print('|' in link_split[0])
                        # print("warning - multiple link options in left side link")

                        if primary_id not in rejected_structures:
                            rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                        rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                        rejected_structures[primary_id][ALTERNATIVE_LIP_IDX] = True
                    # Reject cases where the are multiple alternative links (right side)
                    if '|' in link_split[1]:
                        #print('|' in link_split[1])
                        #print("warning - multiple link options in right side link")
                        if primary_id not in rejected_structures:
                            rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                        rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                        rejected_structures[primary_id][ALTERNATIVE_LIP_IDX] = True
                    # Reject cases where there's uncertainty over the carbon number (left side)
                    if '?' in link_split[0]:
                        #print('?' in link_split[0])
                        #print("warning - uncertain carbon number in left side link")

                        if primary_id not in rejected_structures:
                            rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                        rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                        rejected_structures[primary_id][UNCERTAIN_CARBON_NUM_IDX] = True
                    # Reject cases where there's uncertainty over the carbon number (right side)
                    if '?' in link_split[1]:
                        #print('?' in link_split[1])
                        #print("warning - uncertain carbon number in right side link")
                        if primary_id not in rejected_structures:
                            rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                        rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                        rejected_structures[primary_id][UNCERTAIN_CARBON_NUM_IDX] = True
                # Reject single-sided links
                else:
                    single_sided_link_count += 1
                    # print("WARNING - STRANGE LINK FORMAT")
                    # print("link:", link)
                    # print("PrimaryId:", primary_id)
                    # print(result["Sequence"]["value"])

                    if primary_id not in rejected_structures:
                        rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]

                    rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                    rejected_structures[primary_id][SINGLE_SIDED_LINK_IDX] = True

            #print("link_list:", link_list)

        #print(len(wurcs2_unique_res_list))
        #print(num_unique_residues)
        #print(len(wurcs2_unique_res_list) != num_unique_residues)

        if len(wurcs2_unique_res_list) != num_unique_residues:
            #print("num unique residues found not equal to num_unique_residues)")
            sys.exit(0)

        #print("----------")

        # Only proceed if the glycan wasn't rejected above
        if primary_id not in rejected_structures:
            # Final check - reject if the number of residues isn't one more than the number of links
            if len(wurcs2_res_list) != (len(link_list)+1):
                # print(wurcs2_res_list)
                # print(link_list)
                rejected_structures[primary_id] = ['', False, False, False, False, False, False, False]
                rejected_structures[primary_id][WURCS2_STR_IDX] = wurcs_string
                rejected_structures[primary_id][OTHER_LINK_PROBLEM_IDX] = True
            # Otherwise, all checks passed, so create our unique string representation each reaction, and add to bag of
            # reactions if not already present
            else:

                for link in link_list:
                    link_split = link.split('-')
                    child_res = ord(link_split[0][0].lower())-96
                    child_carbon = link_split[0][1:]
                    parent_res = ord(link_split[1][0].lower())-96
                    parent_carbon = link_split[1][1:]

                    # print("wurcs2_res_list:", wurcs2_res_list)
                    # print("child_res:", child_res)
                    # print("parent_res:", parent_res)
                    # print("link_list:", link_list)

                    reaction = wurcs2_res_list[child_res-1]+'£'+child_carbon+'£'+parent_carbon+wurcs2_res_list[parent_res-1]

                    if primary_id not in reactions_bag:
                        reactions_bag[primary_id] = {'reactions': [reaction]}
                    else:
                        # if reaction not in reactions_bag[primary_id]:
                        if reaction not in reactions_bag[primary_id]['reactions']:
                                reactions_bag[primary_id]['reactions'].append(reaction)
                        else:
                            # print("already in the bag")
                            pass

                    # Also add the motifs (as given in glytoucan db) for this glycan to the reactions bag
                    # (these will be used for comparison with the bag of reactions)
                    if "MotifPrimaryId" in result:
                        motif_id = result["MotifPrimaryId"]["value"]
                        if "motifs" not in reactions_bag[primary_id]:
                            reactions_bag[primary_id]["motifs"] = [motif_id]
                        else:
                            if motif_id not in reactions_bag[primary_id]["motifs"]:
                                reactions_bag[primary_id]["motifs"].append(motif_id)
                    else:
                        reactions_bag[primary_id]["motifs"] = ['None']

    # Print some metrics
    print("single_sided_link_count:", single_sided_link_count)
    print("double_sided_link_count:", double_sided_link_count)
    print("total_link_count:", total_link_count)
    print("single_sided_link_count + double_sided_link_count:", single_sided_link_count + double_sided_link_count)

    print("len(reactions_bag):", len(reactions_bag))

    return reactions_bag