from mimetypes import guess_type
from openai import AzureOpenAI
import json
import os
from tqdm import tqdm
import pandas as pd
import base64
from mimetypes import guess_type
from stage1_prompt import *
import glob
from image_augment import augment_images
import xml.etree.ElementTree as ET
from PIL import Image
import re 
import queue 
import time
import random
import json
import concurrent
top_queue=queue.Queue()
side_normal_queue = queue.Queue()
side_anomaly_queue=queue.Queue()
mix_queue=queue.Queue()
transfer_queue=queue.Queue()

queue_dc={"top":top_queue,"side_normal":side_normal_queue,
          "side_anomaly":side_anomaly_queue,"mix":mix_queue,
          "transfer":transfer_queue
          }
# Function to encode a local image into data URL（GPT4o inference data type） 
def local_image_to_data_url(image_path):
    # Guess the MIME type of the image based on the file extension
    mime_type, _ = guess_type(image_path)
    if mime_type is None:
        mime_type = 'application/octet-stream'  # Default MIME type if none is found

    # Read and encode the image file
    with open(image_path, "rb") as image_file:
        base64_encoded_data = base64.b64encode(image_file.read()).decode('utf-8')

    # Construct the data URL
    return f"data:{mime_type};base64,{base64_encoded_data}"


#Similar to 3D-LLM instruction set construction, the user needs to input the OpenAI API and configure the engine, etc.
def GPT4o(system_prompt, instruction, image):
    data_url = local_image_to_data_url(image)
    messages = [
        { "role": "system", "content": system_prompt },
        { "role": "user", "content": [  
            { 
                "type": "text",
                "text": instruction 
            },
            { 
                "type": "image_url",
                "image_url": {
                    "url": data_url
                }
            }
        ] } 
    ]
    
    response = client.chat.completions.create(
        model="gpt-4o", 
        response_format={"type": "json_object"},
        messages=messages,
        temperature=1.7,
        top_p=0.95,
        max_tokens=16384
    )
    # 生成结果
    output = (response.choices[0].message.content)
    return output

def check_or_initialize_json(json_filename):
    """
    This function checks whether the specified JSON file exists. If not, it initializes the file with empty data for two stages. 
    If the file already exists, it loads the content and returns it.

    Parameters:
    json_filename (str): The name of the JSON file to check or initialize.

    Returns:
    dict: The content of the loaded or initialized JSON file.
    """
    if not os.path.exists(json_filename):
        # If the JSON file does not exist, initialize a new JSON file
        stage1_data_json = {"single":{"Description":[],"QA":[]},
                    "multi":{"Description":[],"QA":[]},
                    "exception":{"Description":[],"QA":[]},
                    "exist_tip":{"Description":[],"QA":[]},
                    "no_tip":{"Description":[],"QA":[]},
                    "50":{"Description":[],"QA":[]},
                    "200":{"Description":[],"QA":[]},
                    "open":{"Description":[],"QA":[]},
                    "close":{"Description":[],"QA":[]},
                    "opening_or_closing":{"Description":[],"QA":[]},                           
                    }
        
        stage2_data_json = {
            "aspirate":{"QA":[]},
            "layout_analysis":{"QA":[]},
            "tip_damaged":{"QA":[]},
            "pipette_type":{"QA":[]},
            "tip_type":{"QA":[]},
            "droptip":{"QA":[]},
            "pickuptip":{"QA":[]},
            "thermocycler_close" :{"QA":[]},
            "thermocycler_open":{"QA":[]}
        }

        with open(json_filename, 'w', encoding='utf-8') as json_file:
            json.dump(stage1_data_json, json_file, ensure_ascii=False, indent=4)
        print(f"{json_filename} does not exist. The JSON file has been initialized.")

        with open(json_filename, 'w', encoding='utf-8') as json_file:
            json.dump(stage2_data_json, json_file, ensure_ascii=False, indent=4)
        print(f"{json_filename} does not exist. The JSON file has been initialized.")

    else:
        # If the JSON file exists, load its content
        with open(json_filename, 'r', encoding='utf-8') as json_file:
            data_json = json.load(json_file)
        print(f"{json_filename} already exists. The content is as follows:")
        print(json.dumps(data_json, ensure_ascii=False, indent=4))

        # If the JSON file exists, load its content
        with open(json_filename, 'r', encoding='utf-8') as json_file:
            data_json = json.load(json_file)
        print(f"{json_filename} already exists. The content is as follows:")
        print(json.dumps(data_json, ensure_ascii=False, indent=4))

    return data_json

def add_to_json(data, json_file_name, module_name, field, content):
    if module_name in data and field in data[module_name]:
        data[module_name][field].append(content)
        with open(json_file_name, 'w', encoding='utf-8') as json_file:
            json.dump(data, json_file, ensure_ascii=False, indent=4)
        print(f"Updated {field} field of {module_name}.")
    else:
        print(f"Error: Invalid module name or field.")

def get_crop_data(crop_data_path):
    crop_dir_list = glob.glob(crop_data_path+"/*")
    for crop_dir in crop_dir_list:
        print(crop_dir)
        crop_file_list = glob.glob(crop_dir+'/crop/*')
        print(crop_file_list)
         


def crop_class_data(annotions_xml,output_csv_path,crop_save_path, img_path, type_):
    """
    Parameters:
    annotions_xml (str): The path to the XML file containing image annotations.
    output_csv_path (str): The path to the CSV file where cropped image information will be saved.
    crop_save_path (str): The folder path where cropped images will be stored.
    img_path (str): The folder path where the original images are stored.
    type_ (str): A string that specifies the image type, either 'side' or other types.

    Returns:
    None

    This function crops images based on the provided annotation XML file and saves the cropped images. It also updates a CSV file with the paths and classes of the cropped images. The function groups and filters images by their class labels, crops the specified regions, and saves them. If the cropped image already exists, the function skips that image.
    """
    if type_=='side':
        type_dict = {"multi":'0','single':'0',
                'exist_tip':'2','200':'2',
                '50':'2','damaged':'2','grasp_anomaly':'2',
                'aspiration_anomaly':'2','aspiration_normally':'2','no_tip':'0'}

    else:
        type_dict = {"1":'1','2':'1',
                '3':'1','4':'1',
                '5':'1','6':'1','7':'1',
                '8':'1','9':'3',
                '10':'3','11':'1','12':'2','trash':'2','close':'3','open':'3','opening_or_closing':'3',"empty":"2"}


    if os.path.exists(output_csv_path):
        df = pd.read_csv(output_csv_path)
    else:
        df = pd.DataFrame(columns=['raw_img_path', 'crop_img_path', 'img_class'])
    # Load the existing crop_img_path into a set to speed up lookups
    existing_crops = set(df['crop_img_path'].values)

    if not os.path.exists(crop_save_path):
        os.makedirs(crop_save_path)

    anomaly_annotations_tree = ET.parse(annotions_xml)
    anomaly_root = anomaly_annotations_tree.getroot()

    for image in anomaly_root.findall('image'):
        if type_=='side':
            class_dict={'0':[],'1':[],'2':[]}
        else:
            class_dict={'0':[],'1':[],'2':[],'3':[]}

        if image.findall('box')==[]:
            continue
        
        new_class_dict={'0':[],'1':[],'2':[],'3':[]}
        for box in image.findall('box'):
            label = box.get('label')
            side_type = str(type_dict[label])
            image_name = image.get('name')
            xtl = int(float(box.get('xtl')))
            ytl = int(float(box.get('ytl')))
            xbr = int(float(box.get('xbr')))
            ybr = int(float(box.get('ybr')))
            class_dict[side_type].append((label, img_path+'/'+image_name,(xtl,ytl,xbr,ybr)))#img_path+'/'+image_name是一样的

        for k,v in class_dict.items():
            if k=='1':
                new_class_dict[k]=v
            elif k=='0':
                new_class_dict[k]+=[x for x in v if x[0] in ['multi','single','no_tip']]
            elif k=='2':
                new_class_dict[k]+=[x for x in v if x[0] in ['exist_tip','200','50','damaged','grasp_anomaly','aspiration_anomaly','aspiration_normally','trash']]
            else:
                new_class_dict[k]+=[x for x in v if x[0] in ['opening_or_closing','open','close']]


        for k,v in new_class_dict.items():
            if v==[]:
                continue
            if k=='2' or k=='0':
                #merging 169_200_exist_tip_grasp_anomaly_damaged中的200_exist_tip_grasp_anomaly_damaged
                cls_name = ','.join([x[0] for x in v])
                image_path = v[0][1]
                coord = v[0][2]
                xtl, ytl, xbr, ybr = coord
                image_path = os.path.join(image_path)
                img = Image.open(image_path)
                cropped_img_path = os.path.join(crop_save_path, f"{image_path.split('/')[-1].split('.')[0]}_{cls_name.replace(',','_')}.jpg")
                # If the cropped image already exists in the CSV, skip this record
                if cropped_img_path in existing_crops:
                    print(f"Cropped image {cropped_img_path} already processed, skipping.")
                    continue
                cropped_img = img.crop((xtl, ytl, xbr, ybr))         
                cropped_img.save(cropped_img_path)
                print(f"Cropped image saved: {cropped_img_path}")
                # Add the new record to the DataFrame
                df = pd.concat([df, pd.DataFrame({
                    'raw_img_path': [image_path],
                    'crop_img_path': [cropped_img_path],
                    'img_class': [cls_name]
                })], ignore_index=True)
                #Update the set of existing cropped images
                existing_crops.add(cropped_img_path)

            else:
                for info_ in v:
                    # 169_multi中的multi
                    cls_name = info_[0]
                    image_path = info_[1]
                    coord = info_[2]
                    xtl, ytl, xbr, ybr = coord
                    #continue
                    image_path = os.path.join(image_path)
                    img = Image.open(image_path)
                    cropped_img_path = os.path.join(crop_save_path, f"{image_path.split('/')[-1].split('.')[0]}_{cls_name.replace(',','_')}.jpg")
                    print(cropped_img_path)
                    print(cls_name)
                    # If the cropped image already exists in the CSV, skip this record
                    if cropped_img_path in existing_crops:
                        print(f"Cropped image {cropped_img_path} already processed, skipping.")
                        continue
                    cropped_img = img.crop((xtl, ytl, xbr, ybr))         
                    cropped_img.save(cropped_img_path)
                    print(f"Cropped image saved: {cropped_img_path}")
                    # Add the new record to the DataFrame
                    df = pd.concat([df, pd.DataFrame({
                        'raw_img_path': [image_path],
                        'crop_img_path': [cropped_img_path],
                        'img_class': [cls_name]
                    })], ignore_index=True)
                    #Update the set of existing cropped images
                    existing_crops.add(cropped_img_path)

    # Save the updated CSV file
    df.to_csv(output_csv_path, index=False)
    print(f"Updated CSV saved to {output_csv_path}")  

def crop_class_data2(annotions_xml,output_csv_path,crop_save_path, img_path, type_):
    """
    Parameters:
    annotions_xml (str): The path to the XML file containing image annotations.
    output_csv_path (str): The path to the CSV file where cropped image information will be saved.
    crop_save_path (str): The folder path where cropped images will be stored.
    img_path (str): The folder path where the original images are stored.
    type_ (str): The image type that determines how the labels are classified.

    Returns:
    None

    This function crops images based on the provided annotation XML file and saves the cropped images. It also updates a CSV file with the paths and classes of the cropped images. The function processes the images based on predefined labels, crops specified regions, and saves them. If the cropped image already exists, the function skips that image.
    """
    dc={"mix正常":"mix_normally",
            "mix异常":"mix_anomaly",
            
            "深孔底部有气泡": "bubble",
            "磁珠不够": "insufficient",
            "磁珠过多": "excessive"
        }
    type_dict ={"mix_normally":'0','mix_anomaly':'0',
                
                'bubble':'0',
                'insufficient':'0','excessive':'0'}
    
    if os.path.exists(output_csv_path):
        df = pd.read_csv(output_csv_path)
    else:
        df = pd.DataFrame(columns=['raw_img_path', 'crop_img_path', 'img_class'])
    # Load the existing crop_img_path into a set to accelerate lookup
    existing_crops = set(df['crop_img_path'].values)

    if not os.path.exists(crop_save_path):
        os.makedirs(crop_save_path)

    anomaly_annotations_tree = ET.parse(annotions_xml)
    anomaly_root = anomaly_annotations_tree.getroot()

    for image in anomaly_root.findall('image'):
        
        class_dict={'0':[]}

        if image.findall('box')==[]:
            continue
        
        new_class_dict={'0':[],'1':[],'2':[],'3':[]}
        for box in image.findall('box'):
            
            image_name = image.get('name')
            xtl = int(float(box.get('xtl')))
            ytl = int(float(box.get('ytl')))
            xbr = int(float(box.get('xbr')))
            ybr = int(float(box.get('ybr')))
            label = box.get('label')

            if label in dc.keys():

                label=dc[label] 

            else:
                label=""
            for tag in image.findall('tag'): 
                if tag.get("label") in dc.keys():
                    temp=","+dc[tag.get("label")]
                    label+=temp
            if label !="" :
                class_dict['0'].append((label, img_path+'/'+image_name,(xtl,ytl,xbr,ybr)))#img_path+'/'+image_name是一样的
 
        for k,v in class_dict.items():
            if v==[]:
                continue
             
            #拼接为169_200_exist_tip_grasp_anomaly_damaged中的200_exist_tip_grasp_anomaly_damaged
            cls_name = v[0][0]
            image_path = v[0][1]
            coord = v[0][2]
            xtl, ytl, xbr, ybr = coord
            image_path = os.path.join(image_path)
            img = Image.open(image_path)
            cropped_img_path = os.path.join(crop_save_path, f"{image_path.split('/')[-1].split('.')[0]}_{cls_name.replace(',','_')}.jpg")
            # If the cropped image already exists in the CSV, skip this record
            if cropped_img_path in existing_crops:
                print(f"Cropped image {cropped_img_path} already processed, skipping.")
                continue
            cropped_img = img.crop((xtl, ytl, xbr, ybr))         
            cropped_img.save(cropped_img_path)
            print(f"Cropped image saved: {cropped_img_path}")
            # Add the new record to the DataFrame
            df = pd.concat([df, pd.DataFrame({
                'raw_img_path': [image_path],
                'crop_img_path': [cropped_img_path],
                'img_class': [cls_name]
            })], ignore_index=True)
            #Update the set of existing cropped images
            existing_crops.add(cropped_img_path) 

    # Save the updated CSV file
    df.to_csv(output_csv_path, index=False)
    print(f"Updated CSV saved to {output_csv_path}")  


def crop(protocol_list):
    """
    Parameters:
    protocol_list (list): The list of protocols, where each protocol contains a path representing different experiments or data sources.

    Returns:
    tuple: A tuple containing three lists: paths to the anomaly crop class CSV files, the normal crop class CSV files, and the top crop class CSV files.

    This function iterates through the protocol list and performs different image cropping tasks depending on the protocol type. The cropping tasks are categorized into side and top normal/anomaly classifications. For bottom protocols, it handles mixed anomaly, normal, and fake classifications. The cropping operations are performed using the `crop_class_data` and `crop_class_data2` functions, and the results are saved to corresponding CSV files.
    """
    anomaly_crop_csv_list=[]
    normally_crop_csv_list=[]
    top_crop_csv_list=[] 
    for single_protocol in protocol_list:
        if "bottom" not in single_protocol:
            protocol_name = single_protocol.split('/')[-1] 
            
            #anomaly-side
            anomaly_csv_path = single_protocol+rf"/side/crop/anomaly_crop_class.csv"
            crop_save_path = single_protocol+rf"/side/crop/anomaly"
            annotions_xml =  single_protocol+rf"/side/raw/{protocol_name}_side_anomaly_annotations.xml"
            anomaly_img_path = single_protocol+rf"/side/raw/anomaly"
            crop_class_data(annotions_xml,anomaly_csv_path,crop_save_path, anomaly_img_path,'side')
            anomaly_crop_csv_list.append(anomaly_csv_path)

            #normally-side
            normally_csv_path = single_protocol+rf"/side/crop/normally_crop_class.csv"
            normally_annotations_xml = single_protocol+rf"/side/raw/{protocol_name}_side_normally_annotations.xml"
            normally_img_path = single_protocol+rf"/side/raw/normally"
            crop_save_path = single_protocol+rf"/side/crop/normally"
            crop_class_data(normally_annotations_xml, normally_csv_path, crop_save_path, normally_img_path,'side')
            normally_crop_csv_list.append(normally_csv_path)

            #top
            top_csv_path = single_protocol+rf"/top/crop/normally_crop_class.csv"
            top_annotations_xml = single_protocol+rf"/top/{protocol_name}_top_annotations.xml"
            top_img_path = single_protocol+rf"/top/raw"
            crop_save_path = single_protocol+rf"/top/crop"
            crop_class_data(top_annotations_xml, top_csv_path, crop_save_path, top_img_path,'top')
            top_crop_csv_list.append(top_csv_path)
        else:
            #mix_anomaly
            anomaly_csv_path = single_protocol+rf"/mix/crop/anomaly_crop_class.csv"
            crop_save_path = single_protocol+rf"/mix/crop/anomaly"
            annotions_xml =  single_protocol+rf"/mix/raw/mix_anomaly_annotations.xml"
            anomaly_img_path = single_protocol+rf"/mix/raw/mix_anomaly"
            crop_class_data2(annotions_xml,anomaly_csv_path,crop_save_path, anomaly_img_path,'bottom')
            anomaly_crop_csv_list.append(anomaly_csv_path)

            #mix_normal
            normally_csv_path = single_protocol+rf"/mix/crop/normal_crop_class.csv"
            crop_save_path = single_protocol+rf"/mix/crop/normally"
            normally_annotations_xml = single_protocol+rf"/mix/raw/mix_normal_annotations.xml"
            normally_img_path = single_protocol+rf"/mix/raw/mix_normal"
            crop_class_data2(normally_annotations_xml, normally_csv_path, crop_save_path, normally_img_path,'bottom')
            normally_crop_csv_list.append(normally_csv_path)

            #mix_fake
            fake_csv_path = single_protocol+rf"/mix/crop/fake_crop_class.csv"
            crop_save_path = single_protocol+rf"/mix/crop/fake"
            fake_annotations_xml = single_protocol+rf"/mix/raw/mix_fake_annotations.xml"
            fake_img_path = single_protocol+rf"/mix/raw/mix_fake"
            crop_class_data2(fake_annotations_xml, fake_csv_path, crop_save_path, fake_img_path,'bottom')
            anomaly_crop_csv_list.append(fake_csv_path)
            
            
            #transfer_anomaly
            anomaly_csv_path = single_protocol+rf"/transfer/crop/anomaly_crop_class.csv"
            crop_save_path = single_protocol+rf"/transfer/crop/anomaly"
            annotions_xml =  single_protocol+rf"/transfer/raw/transfer_anomaly_annotations.xml"
            anomaly_img_path = single_protocol+rf"/transfer/raw/transfer_anomaly"
            crop_class_data2(annotions_xml,anomaly_csv_path,crop_save_path, anomaly_img_path,'bottom')
            anomaly_crop_csv_list.append(anomaly_csv_path)

            #transfer_normal
            normally_csv_path = single_protocol+rf"/transfer/crop/normal_crop_class.csv"
            crop_save_path = single_protocol+rf"/transfer/crop/normally"
            normally_annotations_xml = single_protocol+rf"/transfer/raw/transfer_normal_annotations.xml"
            normally_img_path = single_protocol+rf"/transfer/raw/transfer_normal"
            crop_class_data2(normally_annotations_xml, normally_csv_path, crop_save_path, normally_img_path,'bottom')
            normally_crop_csv_list.append(normally_csv_path)


    return anomaly_crop_csv_list, normally_crop_csv_list, top_crop_csv_list


def augment(csv_data_list):
    
    for csv_data in csv_data_list:
        # Initialize new_map_df with the required columns
        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class'])
        
        # Read the raw mapping DataFrame from the CSV file
        raw_map_df = pd.read_csv(csv_data)

        # Construct the output path and class name
        df_output_path = '/'.join(csv_data.split('/')[:-2])
        class_na = csv_data.split('/')[-1].split('_')[0]
        df_output = df_output_path + f'_{class_na}_augment_map_df.csv'

        # Augment images and update new_map_df
        new_map_df = augment_images(raw_map_df, new_map_df)

        # Check if the output file already exists
        new_map_df.to_csv(df_output, index=False, encoding='utf-8-sig')


def parse_and_save_to_excel(dir_path):
    """
    The function reads the text file from the specified path, uses regular expressions to match configuration parameters, and stores the results in a DataFrame. Each row in the DataFrame represents a configuration and its associated image path. If a specific configuration is found (e.g., "thermocyclerModuleV2"), additional logic is applied to set the appropriate values in the DataFrame.

    Parameters:
    dir_path (str): The path to the text file that contains configuration data and image paths.

    Returns:
    DataFrame: A DataFrame with 12 rows and 2 columns ("配置" and "img_path"), representing each configuration and its corresponding image path.
    """
    flag="None"
    # 初始化DataFrame
    df = pd.DataFrame(index=range(12), columns=['配置',"img_path"])
    df['配置'] = None

    with open(dir_path,"r", encoding='utf-8') as f:
        text=f.read()
    pattern = r'load.*,\s*(\d+)\)$'
    
    ls=[]
    for i in range(1,13):
        ls.append(f"{i}.jpg")
    # 匹配每一行
    lines = text.split('\n')
    lss=[]
    for line in lines:
        match = re.search(pattern, line)
        if match:
            lss.append(line)
    for i in lss:
        pattern = r'["\']([^"]*)["\']'
        match = re.search(pattern, i)
        if match:
            param1 = match.group(1)
            pattern = r',([^)]*)\)'
            match = re.search(pattern, i)
            if match:
                index = int(match.group(1))
                # print(index)
                # 检查索引范围
                if 1 <= index <= 12:
                    if (index==9 and param1=="thermocyclerModuleV2") or (index==10 and param1=="thermocyclerModuleV2"):
                        df.at[9, '配置'] = 'thermocyclerModuleV2'
                        df.at[9,"img_path"]="9-10.jpg"#ls[index-1]
                        df.at[8, '配置'] = 'thermocyclerModuleV2'
                        df.at[8,"img_path"]="9-10.jpg"
                        df.at[11,'配置'] = "trash"
                        df.at[11,'img_path'] = "12.jpg"
                        flag="thermocyclerModuleV2 exist"
                    df.at[index-1, '配置'] = param1
                    df.at[index-1,"img_path"]=ls[index-1]
                    df.at[11,'配置'] = "trash"
                    df.at[11,'img_path'] = "12.jpg"
                    if flag=="thermocyclerModuleV2 exist":
                        df.at[9, '配置'] = 'thermocyclerModuleV2'
                        df.at[9,"img_path"]="9-10.jpg"#ls[index-1]
                        df.at[8, '配置'] = 'thermocyclerModuleV2'
                        df.at[8,"img_path"]="9-10.jpg"
    return df

def bottom_prompt(rows,description_dict,status_requirements_template,status_requirements_template_QA,system_prompt,output_path,tt):
    """
    Function Description:
    Generates experimental phase descriptions and QA results for bottom view images, and saves the results to a global queue.

    Parameters:
    rows (dict): A dictionary containing image paths and classification information, including:
        - 'raw_img_path' (str): Path to the raw image.
        - 'crop_image_path' (str): Path to the cropped image.
        - 'augment_raw_img_path' (str): Path to the augmented raw image.
        - 'augment_crop_img_path' (str): Path to the augmented cropped image.
        - 'img_class' (str): Image classification information, separated by commas.

    description_dict (dict): A dictionary mapping image classes to their corresponding descriptions, used to generate experimental status descriptions.
    status_requirements_template (str): Template for describing the experimental phase.
    status_requirements_template_QA (str): Template for generating QA questions.
    system_prompt (str): System prompt for guiding GPT to generate responses.
    output_path (str): Output path (currently unused, reserved for future extensions).
    tt (str): Task type identifier used to select the corresponding queue.

    Returns:
    None: The results are stored in the global queue `queue_dc`, and no value is returned.

    Functionality:
    1. Extracts image classification information and generates a list of experimental status descriptions.
    2. Generates detailed text for the bottom view experimental status based on the classification descriptions.
    3. Invokes the GPT model to generate descriptive responses and QA answers.
    4. Formats the generated results into a DataFrame and stores them in the specified queue.

    Exception Handling:
    If the model invocation fails 10 consecutive times, prints an error message and stops retrying.
    """

    global queue_dc
    pipette_status=""
    raw_img_path = rows['raw_img_path']
    crop_image_path = rows['crop_image_path']
    augment_raw_img_path =  rows['augment_raw_img_path']
    augment_crop_img_path = rows['augment_crop_img_path']
    img_class = rows['img_class']
    if img_class.startswith(","):
        img_class_list=img_class[1:].split(',')
    else:
        img_class_list = img_class.split(',')
    img_class_description_list = [ description_dict[cl] for cl in img_class_list]
    
    
    #Generate the supplementary bottom view status description based on the class
    for idx, (_, img_description) in enumerate(zip(img_class_list, img_class_description_list)):
        pipette_status+=f"({idx+1}). {img_description}\n"
    print(pipette_status)
 
    horizontal_side_instruction_prompt = f"""This image displays the bottom view of a 96-well biorad_96_wellplate_200ul ,
      where the experiment is conducted using a Magnetic Module V2. The bottom view allows us to determine the current phase of the experiment based 
      on the color of the wells and the height of the liquid within them, indicating whether it is the mix phase or the transfer phase. During the mix phase, it ensures all reaction components are evenly 
        distributed in the solution; for instance, in a PCR reaction, template DNA, primers, dNTPs, buffer, and polymerase must be thoroughly 
        mixed to ensure consistent component ratios in each microliter of reaction mixture. The transfer phase can be used to remove unwanted
          substances or impurities, such as purifying target nucleic acids through multiple washes and transfers during nucleic acid extraction. 
          According to the image, the current experimental phase is the {tt} phase."""
    description_prompt = horizontal_side_instruction_prompt+'\n'+status_requirements_template
    QA_input_prompt = horizontal_side_instruction_prompt+'\n'+status_requirements_template_QA

    horizontal_side_instruction_prompt2 = f"""This image provides a bottom-up view of one row of eight wells from a biorad_96_wellplate_200ul.This row of wells has the following "status": {pipette_status}"""
    description_prompt2 = horizontal_side_instruction_prompt2+'\n'+status_requirements_template
    QA_input_prompt2 = horizontal_side_instruction_prompt2+'\n'+status_requirements_template_QA
    for attempt in range(10):  
        print(f"coming {attempt}")
        second=random.uniform(2, 15)
        time.sleep(second)
        try:
            description_raw_response = GPT4o(system_prompt,description_prompt,raw_img_path)
            qa_raw_response = GPT4o(system_prompt,QA_input_prompt,raw_img_path)

            description_crop_response = GPT4o(system_prompt,description_prompt2,crop_image_path)
            qa_crop_response = GPT4o(system_prompt,QA_input_prompt2,crop_image_path)

            description_raw_aug_response = description_raw_response
            qa_raw_aug_response = qa_raw_response

            description_crop_aug_response = description_crop_response
            qa_crop_aug_response = qa_crop_response

            new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': [img_class],
                            'raw_description':[description_raw_response],'raw_qa':[qa_raw_response],'raw_aug_description':[description_raw_aug_response],'raw_aug_qa':[qa_raw_aug_response],
                            'crop_description':[description_crop_response],'crop_qa':[qa_crop_response],'crop_aug_description':[description_crop_aug_response],'crop_aug_qa':[qa_crop_aug_response]})
            
            queue_dc[tt].put(new_row)
            print("success")
            break
        except Exception as e:
            
            if attempt == 9:  
                print("All attempts failed. Please check the inputs and try again.")

def side_prompt(rows,description_dict,status_requirements_template,status_requirements_template_QA,system_prompt,tt):
    """
    Generate experimental status descriptions and QA results for side view images, 
    and save the results to a global queue.

    Parameters:
    rows (dict): A dictionary containing image paths and classification information, including:
        - 'raw_img_path' (str): Path to the raw image.
        - 'crop_image_path' (str): Path to the cropped image.
        - 'augment_raw_img_path' (str): Path to the augmented raw image.
        - 'augment_crop_img_path' (str): Path to the augmented cropped image.
        - 'img_class' (str): Image classification information, separated by commas.

    description_dict (dict): A dictionary mapping image classes to their corresponding descriptions, used to generate experimental status descriptions.
    status_requirements_template (str): Template for describing the experimental status.
    status_requirements_template_QA (str): Template for generating QA questions.
    system_prompt (str): System prompt for guiding GPT to generate responses.
    tt (str): Task type identifier used to select the corresponding queue.

    Returns:
    None: The results are stored in the global queue `queue_dc`, and no value is returned.

    Functionality:
    1. Extracts image classification information and generates a list of experimental status descriptions.
    2. Generates detailed text for the side view experimental status based on the classification descriptions.
    3. Invokes the GPT model to generate descriptive responses and QA answers.
    4. Formats the generated results into a DataFrame and stores them in the specified queue.

    Exception Handling:
    If the model invocation fails 10 consecutive times, prints an error message and stops retrying.
    """

    global queue_dc
    pipette_status=""
    raw_img_path = rows['raw_img_path']
    crop_image_path = rows['crop_image_path']
    augment_raw_img_path =  rows['augment_raw_img_path']
    augment_crop_img_path = rows['augment_crop_img_path']
    img_class = rows['img_class']
    img_class_list = img_class.split(',')
    img_class_description_list = [ description_dict[cl] for cl in img_class_list]
    
    #Generate the supplementary side view status description based on the class
    for idx, (_, img_description) in enumerate(zip(img_class_list, img_class_description_list)):
        pipette_status+=f"({idx+1}). {img_description}\n"
    print(pipette_status)
    Pipette_module =  "multi-channel(8 channels)"
    horizontal_side_instruction_prompt = f"""The picture shows the Horizontal side view of the alphatool pipette arm, the pipette arm specification is {Pipette_module}. 
    The pipette arm has the following "status": {pipette_status}"""
    description_prompt = horizontal_side_instruction_prompt+'\n'+status_requirements_template
    QA_input_prompt = horizontal_side_instruction_prompt+'\n'+status_requirements_template_QA
    for attempt in range(10):  
        print(f"coming {attempt}")
        second=random.uniform(2, 15)
        time.sleep(second)
        try:
            description_raw_response = GPT4o(system_prompt,description_prompt,raw_img_path)
            qa_raw_response = GPT4o(system_prompt,QA_input_prompt,raw_img_path)

            description_crop_response = GPT4o(system_prompt,description_prompt,crop_image_path)
            qa_crop_response = GPT4o(system_prompt,QA_input_prompt,crop_image_path)

            description_raw_aug_response = description_raw_response
            qa_raw_aug_response = qa_raw_response

            description_crop_aug_response = description_crop_response
            qa_crop_aug_response = qa_crop_response

            new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': [img_class],
                            'raw_description':[description_raw_response],'raw_qa':[qa_raw_response],'raw_aug_description':[description_raw_aug_response],'raw_aug_qa':[qa_raw_aug_response],
                            'crop_description':[description_crop_response],'crop_qa':[qa_crop_response],'crop_aug_description':[description_crop_aug_response],'crop_aug_qa':[qa_crop_aug_response]})
            
            queue_dc[tt].put(new_row)
            print("success")
            break
        except Exception as e:
            
            if attempt == 9:  
                print("All attempts failed. Please check the inputs and try again.")


def get_py_file(path):
    # Get the grandparent directory
    upper_upper_upper_dir = os.path.dirname(os.path.dirname(os.path.dirname(path)))
    
    # Build the target path pattern.
    great_grandparent_dir_name = os.path.basename(upper_upper_upper_dir)
    
    
    return os.path.join(upper_upper_upper_dir,great_grandparent_dir_name+".py")
# protocol_type is cameras
def create_data(protocol_type, df, description_dict, output_path):
    """
    Create descriptions and QA results for experimental images and save them as a CSV file.

    Parameters:
    protocol_type (str): The protocol type, which can be "side", "top", or "bottom", specifying the type of data being processed.
    df (DataFrame): A table containing image paths, cropped paths, augmented paths, and classification information.
    description_dict (dict): A dictionary mapping classification information to descriptive text.
    output_path (str): The path to save the resulting CSV file.

    Returns:
    DataFrame: A table containing image paths, descriptions, and QA results.

    Functionality:
    1. Processes experimental image data based on the protocol type ("side", "top", or "bottom").
    2. Generates descriptions and QA pairs for each image, describing the state of experimental equipment or instruments.
    3. Calls different processing functions (`side_prompt` or `bottom_prompt`) based on the data type.
    4. Uses concurrent processing to accelerate batch data processing.
    5. Saves the generated descriptions and QA results to the specified CSV file.

    Exception Handling:
    - Captures and prints errors in concurrent tasks.
    - Raises `ValueError` if the provided `protocol_type` is invalid.
    system_prompt = "You are an AI assistant that helps people find information."
    """
    description_dict = {str(k): v for k, v in description_dict.items()}
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
    # side
    if protocol_type == "side":
        if "normal" in output_path:
            tt="side_normal"
        else: 
            tt="side_anomaly"
        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class',
                                            'raw_description','raw_qa','raw_aug_description','raw_aug_qa',
                                            'crop_description','crop_qa','crop_aug_description','crop_aug_qa',])
        #description
        status_requirements_template = """Please generate a complete long description of the images based on the above background description and the content seen in the images.The requirements are:
1) The description must have the image seen and cannot be separated from the status;
2) The description should focus on the equipment or instrument itself, and do not speculate on the subsequent experimental content and shooting style.
The json format is as follows(Note the use of double quotes.): {"descriptions":""}"""
        
        #QA
        status_requirements_template_QA="""Please generate two sets of Q&A pairs based on the above background description and the content seen in the images. The requirements are:
1) The questions and answers must be directly related to the images and the answers cannot be separated from the "status";
2) The question can be broad, but the answer should be complete and describe the content of the picture as much as possible.
The json format is as follows (Note the use of double quotes.): {"qa_pairs":{"question":"***","answer":"***"}}"""
        
        max_workers=5
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for idx, row in tqdm(df.iterrows(), total=len(df), desc="Submitting tasks to ThreadPoolExecutor"):
                futures.append(executor.submit(side_prompt, row, description_dict, status_requirements_template, status_requirements_template_QA, system_prompt,tt))
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing side rows in df"):
                try:
                    result = future.result()
                except Exception as e:
                    print(f"An error occurred: {e}")
        df_list=[]
        while not queue_dc[tt].empty():
            df_list.append(queue_dc[tt].get())
        if df_list:
            new_map_df = pd.concat([new_map_df]+ df_list , ignore_index=True)
        new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig') 
        return new_map_df

    #TOP VIEWER
    elif protocol_type == "top":
        print('top data create')
        py_path = get_py_file(df.loc[1]['raw_img_path'])
        top_df = parse_and_save_to_excel(py_path)
        #Get the module for each slot.
        top_dict = {
            str(idx + 1): "usascientific_12_reservoir_22ml" if x['配置'] == "mgi_12_reservoir_22ml" else x['配置']
            for idx, x in top_df.iterrows()
        }
        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class',
                                            'raw_description','raw_qa','raw_aug_description','raw_aug_qa',
                                            'crop_description','crop_qa','crop_aug_description','crop_aug_qa',])
        
        new_rows=[]
        for index, row in df.iterrows():
            img_class = str(row['img_class'])
            
            # check img_class is 'close', 'open', 或 'opening_or_closing'
            if img_class in ['close', 'open', 'opening_or_closing']:
                # Create a new row and modify the img_class value.
                new_row1 = row.copy()
                new_row1['img_class'] = "9"
                new_rows.append(new_row1)

                new_row2 = row.copy()
                new_row2['img_class'] = "10"
                new_rows.append(new_row2)
            else:
                # If img_class is not one of the above values, directly add the original row to the list.
                new_rows.append(row)

        # Create a new DataFrame using all the collected rows at once.
        new_df = pd.DataFrame(new_rows, columns=df.columns) 
        df=new_df
        
        for idx, row in tqdm(df.iterrows(),total=len(df)):
            class_name = str(row['img_class'])
            raw_img_path = row['raw_img_path']
            crop_image_path = row['crop_image_path']
            augment_raw_img_path = row['augment_raw_img_path']
            augment_crop_img_path = row['augment_crop_img_path']
            
            #TODO:check null df
            if idx == 0:
                #This is the image of the AlphaTool experimental layout. Various modules or laboratory consumables can be placed in the layout slots. The slot positions in the image are as follows:
                all_slot_description = 'This is a image of the AlphaTool experiment layout. Various modules or labware may be placed in the layout slots.The slot placement in the figure is as follows:\n'
                #This is used only for counting the thermocyclerModuleV2
                i=0
                #Traverse the current protocol's
                for idx, (slot_num, module) in enumerate(top_dict.items()):
                    if str(module)=='None':
                        module='empty_slot'
                    if module=='thermocyclerModuleV2':
                        i+=1
                    if module=='thermocyclerModuleV2' and i==2:
                        continue
                    img_class_list = [row['img_class'] for idx, row in df.iterrows()]
                    
                    description = description_dict[module]
                    if module=='empty_slot':
                        all_slot_description += f"##slot {idx+1}: (No module or labware) {description}\n"
                    else:
                        all_slot_description += f"##slot {idx+1}: {module} module {description}\n"
                    
                if 'close' in img_class_list:
                    all_slot_description += 'The thermocyclerModuleV2 module has its lid closed, and a thermal cycling experiment may be in progress.\n'
                
                elif 'open' in img_class_list:
                    all_slot_description += 'The thermocyclerModuleV2 module has its lid open and contains the biorad_96_wellplate_200ul_pcr labware, in preparation for subsequent thermal cycling experiments.\n'
                
                elif 'empty' in img_class_list:
                    all_slot_description += 'In addition, the current trash can status is empty, and subsequent experiments can be run with confidence.\n'
                

                top_QA_requirements = """Please generate 5 sets of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must be variety, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description;
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English;
4) The question must include questions and answers about the scene in the picture. For example, "What does this picture describe, and what are the specific slots?".
5) When there are modules or materials, the answer must include the slots corresponding to the modules and materials. Be specific about what is placed in what position.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"}]}"""

                top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
3）When there are modules or materials, the answer must include the slots corresponding to the modules and materials. Be specific about what is placed in what position.
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                description_instruction = all_slot_description + top_description_requirements
                #Perform QA and description generation for the augmented images as well


                qa_instruction = all_slot_description + top_QA_requirements
                
                description_raw_response = GPT4o(system_prompt,description_instruction,raw_img_path)
                qa_raw_response = GPT4o(system_prompt,qa_instruction,raw_img_path)
                print(description_raw_response)
                print(qa_raw_response)
                
                #Generate QA based on the requirements for all blocks in the module
                description_raw_aug_response = GPT4o(system_prompt,description_instruction,augment_raw_img_path)
                qa_raw_aug_response = GPT4o(system_prompt,qa_instruction,augment_raw_img_path)
                print(description_raw_aug_response)
                print(qa_raw_aug_response)

                new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[""], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[""],
                            'img_class': ['top_all'],
                            'raw_description':[description_raw_response],'raw_qa':[qa_raw_response],'raw_aug_description':[description_raw_aug_response],'raw_aug_qa':[qa_raw_aug_response],
                            'crop_description':[''],'crop_qa':[''],'crop_aug_description':[''],'crop_aug_qa':['']})
                new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)

            """Below are the prompts generated for crop"""
            if class_name not in  ['open', 'opening_or_closing']: #[trash, close]
                if class_name in top_dict:
                    # Get the corresponding module.
                    label_name = str(top_dict[class_name])
                else:
                    label_name = class_name

                if label_name == 'None':
                    module_description = description_dict['empty_slot']
                    slot_description = f"The currently displayed slot is slot {class_name}." 
                    top_QA_requirements = """Please generate a Q&A pair based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"}]}"""
                    top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                    top_double_QA_instruction_template = f"""{module_description}{slot_description}\n{alphatool_slot_template}\n{top_QA_requirements}"""

                    top_double_description_instruction_template = f"""{module_description}{slot_description}\n{alphatool_slot_template}\n{top_description_requirements}"""
                    
                    double_description_crop_response = GPT4o(system_prompt,top_double_description_instruction_template,crop_image_path)
                    print(double_description_crop_response)
                    double_qa_crop_response = GPT4o(system_prompt,top_double_QA_instruction_template,crop_image_path)
                    print(double_qa_crop_response)

                    double_description_crop_aug_response = GPT4o(system_prompt,top_double_description_instruction_template,augment_crop_img_path)
                    double_qa_crop_aug_response = GPT4o(system_prompt,top_double_QA_instruction_template,augment_crop_img_path)   
                    print(double_description_crop_aug_response)
                    print(double_qa_crop_aug_response)
                    new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': [class_name],
                            'raw_description':[''],'raw_qa':[''],'raw_aug_description':[''],'raw_aug_qa':[''],
                            'crop_description':[double_description_crop_response],'crop_qa':[double_qa_crop_response],'crop_aug_description':[double_description_crop_aug_response],'crop_aug_qa':[double_qa_crop_aug_response]})
                    new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
                    print(new_map_df)                    

                
                elif (label_name != 'thermocyclerModuleV2') and (label_name != 'None'):
                    if class_name in ['trash','close','empty'] :
                        module_name = label_name
                        module_description = description_dict[module_name]
                        top_QA_requirements = """Please generate 1 set of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"}]}"""
                        top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                        top_double_QA_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
{alphatool_slot_template}\n{top_QA_requirements}"""


                        top_double_description_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
{alphatool_slot_template}\n{top_description_requirements}"""


                        double_description_crop_response = GPT4o(system_prompt,top_double_description_instruction_template,crop_image_path)
                        print(double_description_crop_response)
                        double_qa_crop_response = GPT4o(system_prompt,top_double_QA_instruction_template,crop_image_path)
                        print(double_qa_crop_response)


                        double_description_crop_aug_response = GPT4o(system_prompt,top_double_description_instruction_template,augment_crop_img_path)
                        double_qa_crop_aug_response = GPT4o(system_prompt,top_double_QA_instruction_template,augment_crop_img_path)   
                        print(double_description_crop_aug_response)
                        print(double_qa_crop_aug_response)
                        new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                                "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                                'img_class': [class_name],
                                'raw_description':[''],'raw_qa':[''],'raw_aug_description':[''],'raw_aug_qa':[''],
                                'crop_description':[double_description_crop_response],'crop_qa':[double_qa_crop_response],'crop_aug_description':[double_description_crop_aug_response],'crop_aug_qa':[double_qa_crop_aug_response]})
                        new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
                        print(new_map_df)
                    
                    #OTHER MODULE 
                    elif label_name in ["temperatureModuleV2","magneticModuleV2","heaterShakerModuleV1"] :
                        module_name = label_name
                        labware_name = 'biorad_96_wellplate_200ul_pcr'
                        module_description = description_dict[module_name]
                        labware_description = description_dict[labware_name]
                        # generate 3 sets of Q&A
                        top_QA_requirements = """Please generate 3 sets of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"}]}"""
                        #generate description
                        top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                        top_double_QA_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
2)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_QA_requirements}"""


                        top_double_description_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
2)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_description_requirements}"""


                        
                        double_description_crop_response = GPT4o(system_prompt,top_double_description_instruction_template,crop_image_path)
                        print(double_description_crop_response)
                        double_qa_crop_response = GPT4o(system_prompt,top_double_QA_instruction_template,crop_image_path)
                        print(double_qa_crop_response)
                        double_description_crop_aug_response = GPT4o(system_prompt,top_double_description_instruction_template,augment_crop_img_path)
                        double_qa_crop_aug_response = GPT4o(system_prompt,top_double_QA_instruction_template,augment_crop_img_path)   
                        print(double_description_crop_aug_response)
                        print(double_qa_crop_aug_response)
                        new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                                "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                                'img_class': [class_name],
                                'raw_description':[''],'raw_qa':[''],'raw_aug_description':[''],'raw_aug_qa':[''],
                                'crop_description':[double_description_crop_response],'crop_qa':[double_qa_crop_response],'crop_aug_description':[double_description_crop_aug_response],'crop_aug_qa':[double_qa_crop_aug_response]})
                        new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
                        print(new_map_df)
                        
                    # the experimental consumables separately in the slot
                    else:
                        labware_name = 'biorad_96_wellplate_200ul_pcr' 
                        labware_description = description_dict[labware_name]
                        top_QA_requirements = """Please generate 3 sets of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"}]}"""
                        #Generate a description
                        top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                        top_double_QA_instruction_template = f"""The picture shows the alphatool {labware_name} labware. 
The labware is  to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_QA_requirements}"""


                        top_double_description_instruction_template = f"""The picture shows the alphatool {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_description_requirements}"""


                        
                        double_description_crop_response = GPT4o(system_prompt,top_double_description_instruction_template,crop_image_path)
                        print(double_description_crop_response)
                        double_qa_crop_response = GPT4o(system_prompt,top_double_QA_instruction_template,crop_image_path)
                        print(double_qa_crop_response)
                        double_description_crop_aug_response = GPT4o(system_prompt,top_double_description_instruction_template,augment_crop_img_path)
                        double_qa_crop_aug_response = GPT4o(system_prompt,top_double_QA_instruction_template,augment_crop_img_path)   
                        print(double_description_crop_aug_response)
                        print(double_qa_crop_aug_response)
                        new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                                "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                                'img_class': [class_name],
                                'raw_description':[''],'raw_qa':[''],'raw_aug_description':[''],'raw_aug_qa':[''],
                                'crop_description':[double_description_crop_response],'crop_qa':[double_qa_crop_response],'crop_aug_description':[double_description_crop_aug_response],'crop_aug_qa':[double_qa_crop_aug_response]})
                        new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)

                    
                    
            elif label_name == 'thermocyclerModuleV2' :

                module_name = 'thermocyclerModuleV2'
                labware_name = 'biorad_96_wellplate_200ul_pcr'
                module_description = description_dict[module_name]
                labware_description = description_dict[labware_name]
                
                top_QA_requirements = """Please generate 3 sets of Q&A pairs based on the above background description and the content seen in the picture. Requirements:
1) The questions must include a variety of questions, such as what experimental labware or module is in the picture; it can also be a combination of one or more of background knowledge, layout requirements, possible experiments, etc. Note that the experimental equipment and modules are not strongly bound;
2) The question must be directly related to the picture, and the answer cannot be separated from the background knowledge. In addition, the question is strongly related to the module and weakly related to the labware. The answer must give a complete explanation or description.
3) The wording and grammar of the questions and answers can be varied, but they must all be professional English.
The json format is as follows (note the use of double quotes):{"qa_pairs":[{"question":"***","answer":"***"},{"question":"***","answer":"***"}]}"""
                top_description_requirements = """Please generate one complete long description of the images based on the above background description and the content seen in the images. The requirements are: 
1) The description must have the image seen and cannot be separated from the background knowledge; 
2) The description must include complete background knowledge, layout requirements, possible experiments and other key information, etc.; 
The json format is as follows(Note the use of double quotes.): {"description":""}"""

                top_double_QA_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
2)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_QA_requirements}"""
                print(top_double_QA_instruction_template)

                top_double_description_instruction_template = f"""The picture shows the alphatool {module_name} module and {labware_name} labware. 
The labware is stuck on the module to carry the experimental solution, and the module performs the corresponding experimental operations. 
Their background descriptions are as follows:\n
1)The module background of {module_name} is as follows：{module_description}\n
2)The labware background of {labware_name} is as follows：{labware_description}\n{alphatool_slot_template}\n{top_description_requirements}"""
                print(top_double_description_instruction_template)
                
                double_description_crop_response = GPT4o(system_prompt,top_double_description_instruction_template,crop_image_path)
                print(double_description_crop_response)
                double_qa_crop_response = GPT4o(system_prompt,top_double_QA_instruction_template,crop_image_path)
                print(double_qa_crop_response)

                double_description_crop_aug_response = GPT4o(system_prompt,top_double_description_instruction_template,augment_crop_img_path)
                double_qa_crop_aug_response = GPT4o(system_prompt,top_double_QA_instruction_template,augment_crop_img_path)   
                print(double_description_crop_aug_response)
                print(double_qa_crop_aug_response)
                new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                        "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                        'img_class': [class_name],
                        'raw_description':[''],'raw_qa':[''],'raw_aug_description':[''],'raw_aug_qa':[''],
                        'crop_description':[double_description_crop_response],'crop_qa':[double_qa_crop_response],'crop_aug_description':[double_description_crop_aug_response],'crop_aug_qa':[double_qa_crop_aug_response]})
                new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
                
                

            else:
                print(f'{class_name} is error type')
                exit()
        
        new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        print(f'Complete save {output_path}')
        #exit()
        return new_map_df

    elif protocol_type == "bottom": 
        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class',
                                            'raw_description','raw_qa','raw_aug_description','raw_aug_qa',
                                            'crop_description','crop_qa','crop_aug_description','crop_aug_qa',])
        #description
        status_requirements_template = """Please generate a complete long description of the images based on the above background description and the content seen in the images.The requirements are:
1) The description must have the image seen and cannot be separated from the status;
2) The description should focus on the equipment or instrument itself, and do not speculate on the subsequent experimental content and shooting style.
The json format is as follows(Note the use of double quotes.): {"descriptions":""}"""
        
        #QA
        status_requirements_template_QA="""Please generate two sets of Q&A pairs based on the above background description and the content seen in the images. The requirements are:
1) The questions and answers must be directly related to the images and the answers cannot be separated from the "status";
2) The question can be broad, but the answer should be complete and describe the content of the picture as much as possible.
The json format is as follows (Note the use of double quotes.): {"qa_pairs":{"question":"***","answer":"***"}}"""
        max_workers=50
        if "mix" in output_path:
            tt="mix"
        else:
            tt="transfer"
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = []
            for idx, row in tqdm(df.iterrows(), total=len(df), desc="Submitting tasks to ThreadPoolExecutor"):
                futures.append(executor.submit(bottom_prompt, row, description_dict, status_requirements_template, status_requirements_template_QA, system_prompt,output_path,tt))
            
            # Wait for all tasks to complete
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Processing side rows in df"):
                try:
                    result = future.result()
                except Exception as e:
                    print(f"An error occurred: {e}")
        df_list=[]
        while not queue_dc[tt].empty():
            df_list.append(queue_dc[tt].get())
        if df_list:
            new_map_df = pd.concat([new_map_df]+ df_list , ignore_index=True)
        new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig') 
        return new_map_df

    else:
        raise ValueError("only support bottom、side和top")
    #description_response = GPT4o(system_prompt,description_prompt,img_path)
    #qa_response = GPT4o(system_prompt,qa_prompt,img_path)


def process_protocol(protocol_folder,description_dict):
    """
    Process experimental data in the protocol folder, generate descriptions and QA results based on the experiment type, and save them as CSV files.

    Parameters:
    protocol_folder (str): Path to the protocol folder containing CSV files with experimental data.
    description_dict (dict): A dictionary mapping classification information to descriptive text.

    Functionality:
    1. If the folder path does not contain "bottom", process "side" and "top" data:
    - Read side view anomaly data, side view normal data, and top view data from the folder.
    - Call `create_data` to generate descriptions and QA pairs for each data type and save the results as CSV files.
    - Use concurrent processing to accelerate data generation.
    2. If the folder path contains "bottom", process "mix" and "transfer" data:
    - Read and merge mix data (anomaly, fake, and normal) into a single DataFrame.
    - Call `create_data` to generate descriptions and QA pairs for mix data and save the results as a CSV file.
    - Handle transfer data (currently commented out).
    3. Capture and print any exceptions that occur during task execution.

    Returns:
    None: Results are saved as CSV files, and no value is returned.

    Exception Handling:
    - Captures and prints exceptions from each concurrent task.
    """

    if "bottom" not in protocol_folder:
        augment_side_anomaly_csv = protocol_folder + '/side_anomaly_augment_map_df.csv'
        augment_side_normally_csv = protocol_folder + '/side_normally_augment_map_df.csv'
        augment_top_csv = protocol_folder + '/top_normally_augment_map_df.csv'

        side_anomaly_df = pd.read_csv(augment_side_anomaly_csv)
        side_normally_df = pd.read_csv(augment_side_normally_csv)
        top_df = pd.read_csv(augment_top_csv)

        top_anomaly_output_path =  protocol_folder + '/top_raw_augment_qa_description_map_df.csv'
        
        
        side_anomaly_output_path =  protocol_folder + '/side_anomaly_augment_qa_description_map_df.csv'
        

        side_normally_output_path =  protocol_folder + '/side_normally_augment_qa_description_map_df.csv'
        
        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as thread_executor:
            futures = [
                thread_executor.submit(create_data, 'top', top_df, description_dict, top_anomaly_output_path),
                thread_executor.submit(create_data, 'side', side_anomaly_df, description_dict, side_anomaly_output_path),
                thread_executor.submit(create_data, 'side', side_normally_df, description_dict, side_normally_output_path)
            ]
            # Wait for all tasks to complete
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"{protocol_folder}An error occurred: {e}")
    else:
        a=["mix_anomaly_augment_map_df.csv", "mix_fake_augment_map_df.csv"  ,"mix_normal_augment_map_df.csv"]   

        paths = [os.path.join(protocol_folder, file) for file in a[:3]]
        mixes = [pd.read_csv(path) for path in paths] 
        # Merge three data frames
        combined_mix = pd.concat(mixes, ignore_index=True)

        mix_anomaly_output_path =  protocol_folder + '/mix_raw_augment_qa_description_map_df.csv' 

        b=["transfer_anomaly_augment_map_df.csv","transfer_normal_augment_map_df.csv"] 

        paths = [os.path.join(protocol_folder, file) for file in b[:3]]
        transfers = [pd.read_csv(path) for path in paths] 
        # Merge three data frames
        # combined_transfer = pd.concat(transfers, ignore_index=True) 
        # transfer_anomaly_output_path =  protocol_folder + '/transfer_raw_augment_qa_description_map_df.csv' 

        with concurrent.futures.ThreadPoolExecutor(max_workers=3) as thread_executor:
            futures = [
                thread_executor.submit(create_data, 'bottom', combined_mix, description_dict, mix_anomaly_output_path), 
                # thread_executor.submit(create_data, 'bottom', combined_transfer, description_dict, transfer_anomaly_output_path), 
            ]
            # Wait for all tasks to complete
            for future in concurrent.futures.as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    print(f"{protocol_folder}An error occurred: {e}")
    

def create_stage1_data(protocol_dir_list):
    # get seed description
    description_dict = {}
    df_description = pd.read_excel('./labware_description.xlsx',engine="openpyxl")
    for idx, rows in df_description.iterrows():
        print(rows['labware'],rows['description'])
        description_dict[rows['labware']]=rows['description']
    threads = []
    # Create a thread for each protocol folder.
    with concurrent.futures.ProcessPoolExecutor(max_workers=len(protocol_dir_list)) as process_executor:
    # Create a dictionary to store the corresponding protocol_folder for each future
        futures_to_folders = {process_executor.submit(process_protocol, protocol_folder, description_dict): protocol_folder for protocol_folder in protocol_dir_list}
        
        # Wait for all tasks to complete
        for future in concurrent.futures.as_completed(futures_to_folders):
            protocol_folder = futures_to_folders[future]  # Get the corresponding protocol_folder
            try:
                future.result()  # Try to get the result. If an exception is thrown during the function call, it will be re-raised here
            except Exception as e:
                print(f"An error occurred while processing {protocol_folder}: {e}")
 
    print("All protocols processed.")
def split_list(lst, n):
    """Split the list lst into n equal sublists and remove any empty sublists"""
    if n <= 0:
        raise ValueError("n must be greater than 0")
    
    k, m = divmod(len(lst), n)
    sub_lists = [lst[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]
    
    # Remove empty lists
    sub_lists = [sub_list for sub_list in sub_lists if sub_list]
    
    return sub_lists

def create_data2(protocol_type, df, description_dict, output_path):
    """
    Processes experimental image data based on the protocol type (side, top, or bottom),
    generates stage-2 QA information, and saves the results to a CSV file.

    Parameters:
    protocol_type (str): The protocol type, which can be "side", "top", or "bottom".
    df (DataFrame): A table containing experimental image paths, cropped paths, augmented paths, and classification information.
    description_dict (dict): A dictionary mapping classification information to descriptive text (used in some protocol types).
    output_path (str): The path to save the generated CSV file.

    Returns:
    DataFrame: A table containing image paths, descriptions, and QA results.

    Functionality:
    1. Side View:
    - Analyzes side view data to generate the status of the AlphaTool robotic arm and pipette tips.
    - Produces QA pairs related to experimental equipment and saves the results to the specified CSV file.

    2. Top View:
    - Analyzes top view data to evaluate experimental layout and device states (e.g., thermal cycler open/closed status).
    - Generates QA pairs for layout and thermal cycler state, then saves the results as a CSV file.

    3. Bottom View:
    - Analyzes bottom view data to assess solution, mixing, transfer, magnetic bead retention, and bubble states.
    - Generates QA pairs for experimental statuses and saves the results to a CSV file.

    Exception Handling:
    - Prints error messages for unknown classifications.
    - Captures and logs exceptions during result generation.

    File Saving:
    - The generated DataFrame is saved as a CSV file at the specified `output_path` with UTF-8 encoding.
    """

    system_prompt = "You are an AI assistant that helps people find information."

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

    if protocol_type == "side":
        question = "The image shows the side view of the AlphaTool robotic arm. Please analyze and determine the state of the robotic arm and pipette tips in the field of view."

        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class','stage2_qa'])

        for idx, rows in tqdm(df.iterrows(), total = len(df),  desc="Processing side rows in df"):
            
            answer_dict = {"arm_type":["single_channel","multi_channel"],
                           "tip_status":["present", "absent", "none"],
                           "tip_type":["50_ul","200_ul","none"],
                           "aspiration_status":["normal","abnormal","none"],
                           "tip_condition":["intact", "damaged", "none"],
                           "grasp_condition":["normal","abnormal","none"]}
            
            pipette_status=""
            raw_img_path = rows['raw_img_path']
            crop_image_path = rows['crop_image_path']
            augment_raw_img_path =  rows['augment_raw_img_path']
            augment_crop_img_path = rows['augment_crop_img_path']
            img_class = rows['img_class']
            
            img_class_list = img_class.split(',')
            # print('img_class:', img_class_list)
            for class_ in img_class_list:
                #arm_type
                if class_=="single" or class_=="multi":
                    if class_=="single":
                        answer_dict["arm_type"] = "single_channel"
                    else:
                        answer_dict["arm_type"] = "multi_channel"

                # aspiration_status
                elif class_=="aspiration_anomaly" or class_=="aspiration_normally":
                    if class_=="aspiration_anomaly":
                        answer_dict['aspiration_status'] = "abnormal"
                    else:
                        answer_dict['aspiration_status'] = "normal"
                
                # tip_type
                elif class_=='50' or class_=='200' :
                    if class_=='50':
                        answer_dict['tip_type'] = "50_ul"
                    else: # class_=='200'
                        answer_dict['tip_type'] = "200_ul"

                
                # tip_status
                elif class_=='exist_tip' or class_=='no_tip':
                    if class_=='exist_tip':
                        answer_dict['tip_status'] = "present"
                    else: # class_=='200'
                        answer_dict['tip_status'] = "absent"
                
                # tip_condition
                elif class_=='damaged':# or class_=='no_tip':
                    answer_dict['tip_condition'] = "damaged"
                
                elif class_=="grasp_anomaly":
                    answer_dict['grasp_condition'] = "abnormal"
                
                else:
                    print(f"error class: {class_}")
                    
            if 'exist_tip' not in img_class_list:
                answer_dict['tip_status'] = "none"
                answer_dict['tip_type'] = "none"
                answer_dict['aspiration_status'] = "none"
                answer_dict['tip_condition'] = "none"
                answer_dict['grasp_condition'] = "none"
            else:
                if 'grasp_anomaly' not in img_class_list:
                    answer_dict['grasp_condition'] = "normal"
                if 'damaged' not in img_class_list:
                    answer_dict['tip_condition'] = "intact"
            
            if "single" not in img_class_list or "multi" not in img_class_list:
                answer_dict['arm_type'] = "multi_channel"
                
            if "aspiration_anomaly" not in img_class_list or "aspiration_normally" not in img_class_list:
                answer_dict['aspiration_status'] = "none" 
            
            dcb={}
            dcb['Q']=question
            dcb["A"]=answer_dict
            
            stage2_qa=dcb
            new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': [img_class],
                            'stage2_qa':[stage2_qa],})
                            #'raw_description':[description_raw_response],'raw_qa':[qa_raw_response],'raw_aug_description':[description_raw_aug_response],'raw_aug_qa':[qa_raw_aug_response],
                            #'crop_description':[description_crop_response],'crop_qa':[qa_crop_response],'crop_aug_description':[description_crop_aug_response],'crop_aug_qa':[qa_crop_aug_response]})
            new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
            # print(new_map_df)
            new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig')
            print(f'Complete save {output_path}')

        return new_map_df
      

    #TOP VIEWER
    elif protocol_type == "top":
        #open，close，opening_or_closing，empty，trash，empty_slot
        print('top data create')
        print(df)
        py_path = '/'.join([os.getcwd()]+df.loc[1]['raw_img_path'].split('/')[1:3] + [df.loc[1]['raw_img_path'].split('/')[2]+'.py'])
       
        top_df = parse_and_save_to_excel(py_path)
        print(top_df)
        

        top_dict = {str(idx+1):x['配置'] for idx,x in top_df.iterrows()}
        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class','layout_qa','thermocycler_qa'])


        for idx, row in tqdm(df.iterrows(),total=len(df)):
            class_name = str(row['img_class'])
            raw_img_path = row['raw_img_path']
            crop_image_path = row['crop_image_path']
            augment_raw_img_path = row['augment_raw_img_path']
            augment_crop_img_path = row['augment_crop_img_path']

            img_class_list = class_name.split(',')
            
            #all top view
            if idx == 0:
                #print(top_dict)
                layout_question="Please analyze the information for all slot positions of the AlphaTool shown in the image."
                print(top_df)
                layout_answer_dict={}
                
                for idx,row in top_df.iterrows():
                    temp=f'slot {idx+1}'
                    layout_answer_dict[temp]=row['配置']
                
                thermocycler_question = "Please determine the status of the AlphaTool thermal cycler during its operation."
                if 'close' in img_class_list:
                    thermocycler_answer = {"thermocycler_status":"closed"} 
                elif 'open' in img_class_list:
                    thermocycler_answer = {"thermocycler_status":"open"}
                else:
                    thermocycler_answer = {"thermocycler_status":"none"}
                
                aaa=json.dumps({"Q":thermocycler_question,"A":thermocycler_answer})
                
                bbb=json.dumps({"Q":layout_question,"A":layout_answer_dict})
                dc={"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': ['top_all'],'layout_qa':[bbb],
                         'thermocycler_qa':[aaa]    
                            }
                 
                new_row = pd.DataFrame(   dc )
                new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)

            else:
                continue
        
        new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig')
        print(f'Complete save {output_path}')
        return new_map_df

    else:
        question = "The image shows a bottom view of either an 8-well section or the entire 96-well plate. Please analyze the image to determine the status of the wells within the field of view."

        new_map_df = pd.DataFrame(columns=['raw_img_path','crop_image_path', 'augment_raw_img_path', 'augment_crop_img_path', 'img_class','stage2_qa'])
        type_dict ={"mix_normally":'0','mix_anomaly':'0',
                'retention_normal':'0','retention_anomaly':'0',
                'transfer_normal':'0','transfer_anomaly':'0',
                'bubble':'0',
                'insufficient':'0','excessive':'0'}
        for idx, rows in tqdm(df.iterrows(), total = len(df),  desc="Processing side rows in df"):
            answer_dict = {
                           "solution_magnetic_beads_status":["insufficient", "excessive", "none"], 
                           "mix_status":["normal","abnormal","none"],
                           "transfer_status":["normal","abnormal","none"],
                           "magnetic_beads_retention":["normal","abnormal","none"],
                           "bubble":["exist","none"], 
                           }
            ddc         ={
                           "solution_magnetic_beads_status":"none", 
                           "mix_status":"none",
                           "transfer_status":"none",
                           "magnetic_beads_retention":"none",
                           "bubble":"none", 
                           }
            raw_img_path = rows['raw_img_path']
            crop_image_path = rows['crop_image_path']
            augment_raw_img_path =  rows['augment_raw_img_path']
            augment_crop_img_path = rows['augment_crop_img_path']
            img_class = rows['img_class']
            img_class_list = img_class.split(',')
            for cls in img_class_list:
                # 溶液颜色
                if cls in answer_dict["solution_magnetic_beads_status"]:
                    ddc["solution_magnetic_beads_status"]=cls 
                # mix
                if cls =="mix_normally" or cls =="mix_anomaly":
                    if cls =="mix_normally":
                        ddc["mix_status"]="normal"
                    else:
                        ddc["mix_status"]="abnormal" 
                #  transfer
                if cls =="transfer_normal" or cls =="transfer_anomaly":
                    if cls =="transfer_normal":
                        ddc["transfer_status"]="normal"
                    else:
                        ddc["transfer_status"]="abnormal" 
                # 磁珠保留
                if cls =="retention_normal" or cls =="retention_anomaly":
                    if cls =="retention_normal":
                        ddc["magnetic_beads_retention"]="normal"
                    else:
                        ddc["magnetic_beads_retention"]="abnormal" 
                # 底部气泡
                if cls =="bubble":
                   ddc["bubble"]="exist" 
            dcb={}
            dcb['Q']=question
            dcb["A"]=ddc
            stage2_qa=dcb
            new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_image_path], 
                            "augment_raw_img_path":[augment_raw_img_path],"augment_crop_img_path":[augment_crop_img_path],
                            'img_class': [img_class],
                            'stage2_qa':[stage2_qa]})
                            #'raw_description':[description_raw_response],'raw_qa':[qa_raw_response],'raw_aug_description':[description_raw_aug_response],'raw_aug_qa':[qa_raw_aug_response],
                            #'crop_description':[description_crop_response],'crop_qa':[qa_crop_response],'crop_aug_description':[description_crop_aug_response],'crop_aug_qa':[qa_crop_aug_response]})
            new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)
            # print(new_map_df)
            new_map_df.to_csv(output_path, index=False, encoding='utf-8-sig')
            print(f'完成save {output_path}')
 
    #description_response = GPT4o(system_prompt,description_prompt,img_path)
    #qa_response = GPT4o(system_prompt,qa_prompt,img_path)

def process_protocol_folder(protocol_folder, description_dict):
    
    if "bottom" not in protocol_folder:
        augment_side_anomaly_csv = os.path.join(protocol_folder, 'side_anomaly_augment_map_df.csv')
        augment_side_normally_csv = os.path.join(protocol_folder, 'side_normally_augment_map_df.csv')
        augment_top_csv = os.path.join(protocol_folder, 'top_normally_augment_map_df.csv')

        side_anomaly_df = pd.read_csv(augment_side_anomaly_csv)
        side_normally_df = pd.read_csv(augment_side_normally_csv)
        top_df = pd.read_csv(augment_top_csv)

        # Define output paths
        top_anomaly_output_path = os.path.join(protocol_folder, 'stage2_top_raw_augment_qa_description_map_df.csv')
        side_anomaly_output_path = os.path.join(protocol_folder, 'stage2_side_anomaly_augment_qa_description_map_df.csv')
        side_normally_output_path = os.path.join(protocol_folder, 'stage2_side_normally_augment_qa_description_map_df.csv')

        # Sequentially call create_data2 for each dataframe
        create_data2('top', top_df, description_dict, top_anomaly_output_path)
        # print("top")
        create_data2('side', side_anomaly_df, description_dict, side_anomaly_output_path)
        # print("side1") 
        create_data2('side', side_normally_df, description_dict, side_normally_output_path)
        # print("side2")
        print("successfully")
    else:
        mix_csv = os.path.join(protocol_folder, 'mix_raw_augment_qa_description_map_df.csv')
        # transfer_csv = os.path.join(protocol_folder, 'transfer_raw_augment_qa_description_map_df.csv') 

        mix_df = pd.read_csv(mix_csv)
        # transfer_df = pd.read_csv(transfer_csv) 

        # Define output paths 
        mix_output_path = os.path.join(protocol_folder, 'stage2_mix_augment_qa_description_map_df.csv')
        # transfer_output_path = os.path.join(protocol_folder, 'stage2_transfer_augment_qa_description_map_df.csv')

        # Sequentially call create_data2 for each dataframe
        create_data2('mix', mix_df, description_dict, mix_output_path) 
        # create_data2('transfer', transfer_df, description_dict, transfer_output_path)

def create_stage2_data(protocol_dir_list):
    # Get seed description
    df_description = pd.read_excel('./labware_description.xlsx')
    description_dict = {rows['labware']: rows['description'] for _, rows in df_description.iterrows()}

    print(protocol_dir_list)
    from concurrent.futures import ThreadPoolExecutor, as_completed
    # Use ThreadPoolExecutor to run process_protocol_folder in parallel threads
    with ThreadPoolExecutor(max_workers=len(protocol_dir_list)) as executor:  # Adjust max_workers as needed
        futures = {executor.submit(process_protocol_folder, protocol_folder, description_dict): protocol_folder for protocol_folder in protocol_dir_list}
        
        for future in as_completed(futures):
            protocol_folder = futures[future]
            try:
                future.result()  # This will raise an exception if the function did
            except Exception as exc:
                print(f'Processing {protocol_folder} generated an exception: {exc}')

def filp_img(image_path):
    image = Image.open(image_path)
    image = image.rotate(180)
    image.save(image_path)

if __name__=="__main__":
    protocol_path = "./Alphatool_data/"
    protocol_dir_list = [d for d in glob.glob(protocol_path + '/*') if os.path.isdir(d)   ] 
    protocol_dir_list=["./Alphatool_data/bottom/"]
    
    import concurrent.futures
    
    anomaly_crop_csv_list, normally_crop_csv_list, top_crop_csv_list  = crop(protocol_dir_list)  
    
    anomaly_parts = split_list(anomaly_crop_csv_list, 12)
    normal_parts = split_list(normally_crop_csv_list, 12)
    import concurrent
    # Use ThreadPoolExecutor to create multiple threads to process each part
    with concurrent.futures.ThreadPoolExecutor() as executor:
        # For each part of anomaly_crop_csv_list and normally_crop_csv_list a thread is started
        futures = [
            *[
                executor.submit(augment, part)
                for part in anomaly_parts + normal_parts
            ],
            # If top_crop_csv_list also needs to be processed, you can add it here
            executor.submit(augment, top_crop_csv_list)
        ]
        
        # Wait for all tasks to complete
        for future in concurrent.futures.as_completed(futures):
            try:
                result = future.result()  # Get the result or throw an exception
            except Exception as e:
                print(f"An error occurred: {e}") 
    print("crop")
     
    # Only flip the top related graphs
    for file_path in protocol_dir_list:
        top_file_path = os.path.join(file_path,'top' )
        raw_top_file_path =os.path.join(top_file_path,'raw') 
        crop_top_file_path =os.path.join( top_file_path,'crop' )
        augment_path= os.path.join(file_path,'augment' )
        for raw_img in glob.glob(str(raw_top_file_path)+'/*'):
            if 'jpg' in raw_img: 
                filp_img(raw_img) 
        for crop_img  in glob.glob(str(crop_top_file_path)+'/*'):
            if 'jpg' in crop_img:
                filp_img(crop_img) 
        for augment_img  in glob.glob(str(augment_path)+'/*'):
            if 'jpg' in augment_img:
                filp_img(augment_img) 

    #GPT4o server
    #The user is required to provide their own Azure-OpenAI or OpenAI GPT-4 API key and other relevant information.
    client = AzureOpenAI(
        azure_endpoint = "***", 
        api_key="***",  
        api_version="***"
        )
    # Construct QA and descriptive data pools through seed data.
    ## stage1  
    # The two functions create_stage1_data and create_stage2_data are independent. Run whichever function you want and comment out the other one.
    # Note that this is to generate a csv file for QA, which requires connecting to the auzure API. Very slow
    create_stage1_data(protocol_dir_list)

    # This generates the csv of stage2 qa, which is very fast.
    create_stage2_data(protocol_dir_list)