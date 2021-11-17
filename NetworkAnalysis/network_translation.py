import sys, os

"""
    NetworkAnalysis
    2017 Joaquim Aguirre-Plans 
    Structural Bioinformatics Laboratory
    Universitat Pompeu Fabra
"""

def translate(input_network, output_network, translation_file, input_format, output_format, input_nodes, output_nodes):
    """
    Translates the input network (and input_nodes) using the translation file, which contains the translations from
    one type of code to another. Outputs the translated network in a given output_format.
    """


    if not fileExist(input_network):
        print("File with input network is missing\n")
        sys.exit(10)

    if not output_network:
        print("Output network file is not defined!\n")
        sys.exit(10)



    #--------------------------------#
    #   PARSE THE TRANSLATION FILE   #
    #--------------------------------#

    new={}

    # Fill new with the translations
    if fileExist(translation_file):
        ft=open(translation_file,"r")
        for line in ft:
            word=line.strip().split("\t")
            node=word[0]
            new.setdefault(node,set())
            for transet in word[1].split("'"):
                for trans in transet.split(","):
                    if trans != "" and trans != ",":
                        name="_".join([str(x) for x in trans.split()])
                        new[node].add(name)
        ft.close()
        translation=True
    else:
        translation=False



    #---------------------------#
    #   TRANSLATE THE NETWORK   #
    #---------------------------#

    with open(input_network, 'r') as input_network_fd, open(output_network, 'w') as output_network_fd:

        # seeds=set()
        # if fileExist(options.seed):
        #  fn=open(options.seed,'r')
        #  for line in fn:
        #   word=line.split()
        #   seeds.add(word[0])
        #  fn.close()

        nodes=set()

        for line in input_network_fd:

            # Process the input network
            skip=False
            if input_format == 'multi-fields' :
                word=line.strip().split("\t")
                p=word[0] # Node 1
                q=word[1] # Node 2
                try:
                    score=word[2]       # Code for Janet # This is not the score, it is the source, but it is ok
                except:
                    score=1.0
            elif input_format == 'sif':
                word=line.split()
                if len(word)>2:
                    p=word[0]
                    q=word[2]
                    score=1.0
                else:skip=True
            elif input_format == 'raw':
                word=line.split()
                if len(word)>1:
                    p=word[0]
                    q=word[1]
                    score=1.0
                else:skip=True
            else:
                word=line.split()
                p=word[0]
                score=float(word[1])
                q=word[2]
            info="".join([str(word[i]) for i in range(3,len(word))])
            info="\t".join([str(word[i]) for i in range(3,len(word))])       # Code for Janet
            if not translation:
                new.setdefault(p,set()).add(p)
                new.setdefault(q,set()).add(q)

            # Print the output network
            if p in new and q in new:
                for a in new[p]:
                    nodes.add(a)
                    for b in new[q]:
                        nodes.add(b)
                        if output_format == 'multi-fields' :
                            output_network_fd.write("{0}\t{1}\t{2}\t{3}\n".format(a,b,score,info))       # Code for Janet
                        elif output_format == 'sif' :
                            if not skip: output_network_fd.write("{0}\tinteraction\t{1}\n".format(a,b))
                        else:
                            output_network_fd.write("{0} {1:f} {2}\t\t{3:s}\n".format(a,score,b,info))




    #-----------------------------#
    #   TRANSLATE THE NODE FILE   #
    #-----------------------------#

    if fileExist(input_nodes) and fileExist(output_nodes):

        with open(input_nodes, 'r') as input_nodes_fd, open(output_nodes, 'w') as output_nodes_fd:

            if output_format ==  'sif' : output_nodes_fd.write("#%s\n"%(translation_file))

            for line in input_nodes_fd:
                if line[0] == '#':
                    continue
                if input_format == 'multi-fields' :
                    word=line.split()
                    p=word[0]
                    info="\t".join([str(word[i]) for i in range(1,len(word))])
                else:
                    word=line.split()
                    p=word[0]
                    score=float(word[1])
                    if score <= 0.0: score=options.score
                    info=" ".join([str(word[i]) for i in range(2,len(word))])

                # Print the output nodes
                if new.has_key(p):
                    if output_format ==  'sif' :
                        output_nodes_fd.write("{0} = {1}\n".format(p,";".join([x for x in new[p]])))
                    else:
                        for a in new[p]:
                            if output_format == 'multi-fields':
                                output_nodes_fd.write("{0}\t{1}\n".format(a,info))
                            else:
                                # About STRING FORMATTING: 
                                # At the left of ":"  --> 0 is the first format, 1 is the second and 2 is the third
                                # At the right of ":" --> "10" is to align to the left (I have deleted it)
                                #                     --> "f" is to add a float, and it indicates the number of decimal positions after the dot
                                # Good examples at    --> https://pyformat.info/
                                output_nodes_fd.write("{0} {1:0.5f} {2:s}\n".format(a,score,info))

            # else:
            #     for a in nodes:
            #         if a in seeds: score=1.0
            #         else:          score=options.score
            #         if output_format == 'netscore' :
            #             if a is not None: output_nodes_fd.write("{0}\t{1:10.5f}\t{2:10.5f}\t{3:10.5f}\n".format(a,1.,1.,score))
            #         elif output_format ==  'sif' :
            #             if translation:
            #                 if new.has_key(a): 
            #                     if a is not None: output_nodes_fd.write("{0} = {1}\n".format(a,";".join([x for x in new[a]])))
            #                 else:
            #                     if a is not None: output_nodes_fd.write("{0}\n".format(a))
            #             else:
            #                 if a is not None: output_nodes_fd.write("{0}\n".format(a))
            #         else:
            #             if a is not None: output_nodes_fd.write("{0} {1:10.5f}\n".format(a,score))

    return


def fileExist (file):               #Checks if a file exists AND is a file
    if file is not None: return os.path.exists(file) and os.path.isfile(file)
    else: return False
