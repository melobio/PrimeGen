double_instruction_template = "The picture shows the alphatool {module_name} module and {labware_name} labware. The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. Their background descriptions are as follows:\n1)The module background of {module_name} is as follows：{module_description}\n2)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot}\n{requirements}"

instruction_template = "The picture shows the {module_name} module of alphatool, and its background is as follows: {description}\n{alphatool_slot}\n{requirements}"

status_instruction_template = "The picture shows the {module_name} module of alphatool, and its background is as follows: {description}\n{requirements}"
 

requirements_template = """Please generate 5 complete long descriptions of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
3) The description sentences and grammar must be different.
The json format is as follows(Note the use of double quotes.): {"descriptions":[{"description":"***"},{"description":"***"}...]}"""

status_requirements_template = """Please generate 5 complete long descriptions of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge,, possible experiments and other key information, etc.; 
3) The description sentences and grammar must be different.
The json format is as follows(Note the use of double quotes.): {"descriptions":[{"description":"***"},{"description":"***"}...]}
"""
double_requirements_template = """Please generate 5 complete long descriptions of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
3) The description sentences and grammar must be different.
The json format is as follows(Note the use of double quotes.): {"descriptions":[{"description":"***"},{"description":"***"}...]}
"""


QA_requirements_template = """Please generate 10 sets of Q&A pairs based on the above background description and the content seen in the images. The requirements are:
1) The questions must include simple questions, such as what labware or modules are in the images; they can also be about one or a combination of background knowledge, layout requirements, possible experiments, etc.; not every QA needs to cover all aspects completely.
2) The questions and answers must be directly related to the images and the answers cannot be separated from the background knowledge.
3) The wording and grammar of the questions and answers must be different.
The json format is as follows (Note the use of double quotes.): {"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"},...]}"""

status_QA_requirements_template = """Please generate 10 sets of Q&A pairs based on the above background description and the content seen in the images. The requirements are:
1) The questions must include simple questions, such as what labware or modules are in the images; they can also be about one or a combination of background knowledge, possible experiments, etc.; not every QA needs to cover all aspects completely.
2) The questions and answers must be directly related to the images and the answers cannot be separated from the background knowledge.
3) The wording and grammar of the questions and answers must be different.
The json format is as follows (Note the use of double quotes.): {"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"},...]}"""

double_QA_requirements_template = """Please generate 10 sets of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"},...]}"""

alphatool_slot_template = """#Alahtool Slot environment simulation:  
----------------------------
 |  1  |  5  |  9  |
 ---------------------------
 |  2  |  6  |  10  |
 ---------------------------
 |  3  |  7  |  11  |
 ---------------------------
 |  4  |  8  |  12  |
 ---------------------------"""