import os
import random
import json
import pandas as pd
import glob
from functools import reduce

def get_conversation(query, response, img_path, model_type):
    data_type_list_1 = ['glm4v']
    data_type_list_2 = ['phi35', 'MiniCPM','qwen2vl','internvl2', 'llama32']
    data_type_list_3 = ['phi3'] 
    #/home/YJhou/compute/VLM-2024
    img_path=img_path.replace("./Alphatool_data", "/home/Tli/workspace/PCR/Alphatool_data")
    if model_type in data_type_list_1:
        conversation= {"query": query, "response": response, "images": [img_path]}
    elif model_type in data_type_list_2:
        conversation = {"query": f"<image>{query}", "response": response, "images": [img_path]}
    elif model_type in data_type_list_3:
        conversation = {"conversations":[{"from": "user", "value": f"<img>{img_path}</img>{query}"},{"from": "assistant", "value": response}]}
    else:
        print(f"error {model_type} not supported")

    return conversation


def construct_data(data_path,data_type, model_type):
    json_data=[]
    image_set=set()
    global count,count1
    df = pd.read_csv(data_path)
    if data_type == "top":
        for idx,row in df.iterrows():
            try:
                if idx == 0:
                    qa_list = json.loads(row['raw_qa'])['qa_pairs']
                    description = json.loads(row['raw_description'])['description']
                    # 整个版位布局的问答 
                    for qa_ in qa_list:
                        conversation = get_conversation(qa_['question'], qa_['answer'], row['raw_img_path'], model_type)
                        json_data.append(conversation)

                    top_question = 'Please describe the contents of the alphatool experiment layout in the image.'
                    conversation = get_conversation(top_question, description, row['raw_img_path'], model_type)
                    json_data.append(conversation)
                else:

                    crop_description = json.loads(row['crop_description'])['description']
                    crop_qa = json.loads(row['crop_qa'])['qa_pairs'] [0]
                    crop_img = row['crop_image_path']
                    conversation = get_conversation(crop_qa['question'], crop_qa['answer'], crop_img, model_type)
                    json_data.append(conversation)
                    crop_description_query = "Please describe the alphatool device or slot in the image."
                    conversation = get_conversation(crop_description_query, crop_description, crop_img, model_type)
                    json_data.append(conversation)
            except json.decoder.JSONDecodeError:
                count+=1
        return json_data
    
    elif data_type == "side_anomaly":
        for idx,row in df.iterrows():  
            qa_list = json.loads(row['raw_qa'])['qa_pairs'] 
            if(isinstance(qa_list,list)):
                for qa_ in qa_list:
                    conversation = get_conversation(qa_['question'], qa_['answer'], row['raw_img_path'], model_type)
                    json_data.append(conversation) 
            else:
                conversation = get_conversation(qa_list['question'], qa_list['answer'], row['raw_img_path'], model_type)
                json_data.append(conversation) 
            if row['raw_img_path'] not in image_set:
                # print(json.loads(row['raw_description']))
                description = json.loads(row['raw_description'])['descriptions'] 
                side_question = 'Ignore the background content and please provide a detailed description of the Alphatool pipette and its current status in the image. ' 
                conversation = get_conversation(side_question, description, row['raw_img_path'], model_type)
                json_data.append(conversation)
                image_set.add( row['raw_img_path']) 
            crop_description = json.loads(row['crop_description'])['descriptions']
            # print(json.loads(row['crop_qa'])['qa_pairs'])
            crop_img = row['crop_image_path']
            crop_description_query = "Ignore the background content and please provide a detailed description of the thing under the Alphatool pipette and its current status in the image."
            
            conversation = get_conversation(crop_description_query, crop_description, crop_img, model_type)
            json_data.append(conversation)
            crop_qa = json.loads(row['crop_qa'])['qa_pairs'] 
            
            if (isinstance(crop_qa,list)):
                for c in crop_qa:
                    conversation = get_conversation(c['question'], c['answer'], crop_img, model_type)
                    json_data.append(conversation)
            else:
                conversation = get_conversation(crop_qa['question'], crop_qa['answer'], crop_img, model_type)
                json_data.append(conversation)
        return json_data

    elif data_type == "side_normally":
        for idx,row in df.iterrows(): 
            qa_list = json.loads(row['raw_qa'])['qa_pairs'] 
            if(isinstance(qa_list,list)):
                for qa_ in qa_list:
                    conversation = get_conversation(qa_['question'], qa_['answer'], row['raw_img_path'], model_type)
                    json_data.append(conversation) 
            else:
                conversation = get_conversation(qa_list['question'], qa_list['answer'], row['raw_img_path'], model_type)
                json_data.append(conversation) 
            if row['raw_img_path'] not in image_set:
                description = json.loads(row['raw_description'])['descriptions'] 
                side_question = 'Ignore the background content and please provide a detailed description of the Alphatool pipette and its current status in the image. ' 
                conversation = get_conversation(side_question, description, row['raw_img_path'], model_type)
                json_data.append(conversation)
                image_set.add( row['raw_img_path'])

    
            crop_description = json.loads(row['crop_description'])['descriptions']
            crop_img = row['crop_image_path']
            crop_description_query = "Ignore the background content and please provide a detailed description of the thing under the Alphatool pipette and its current status in the image."
            
            conversation = get_conversation(crop_description_query, crop_description, crop_img, model_type)
            json_data.append(conversation)
            crop_qa = json.loads(row['crop_qa'])['qa_pairs'] 
            
            if (isinstance(crop_qa,list)):
                for c in crop_qa:
                    conversation = get_conversation(c['question'], c['answer'], crop_img, model_type)
                    json_data.append(conversation)
            else:
                conversation = get_conversation(crop_qa['question'], crop_qa['answer'], crop_img, model_type)
                json_data.append(conversation)
        return json_data

    else:
        #bottom view
        for idx,row in df.iterrows():
            qa_list = json.loads(row['raw_qa'])['qa_pairs'] 
            if(isinstance(qa_list,list)):
                for qa_ in qa_list:
                    conversation = get_conversation(qa_['question'], qa_['answer'], row['raw_img_path'], model_type)
                    json_data.append(conversation) 
            else:
                conversation = get_conversation(qa_list['question'], qa_list['answer'], row['raw_img_path'], model_type)
                json_data.append(conversation) 
            
            if row['raw_img_path'] not in image_set:
                description = json.loads(row['raw_description'])['descriptions'] 
                side_question = 'Ignore the background content and please provide a detailed description of the 96-well biorad_96_wellplate_200ul and its current status in the image. ' 
                conversation = get_conversation(side_question, description, row['raw_img_path'], model_type)
                json_data.append(conversation)
                image_set.add( row['raw_img_path'])

    
            crop_description = json.loads(row['crop_description'])['descriptions']
            crop_img = row['crop_image_path']
            crop_description_query = "Please describe in detail the current status of the single row of eight wells visible in the image, ignoring the background."
            
            conversation = get_conversation(crop_description_query, crop_description, crop_img, model_type)
            json_data.append(conversation)
            crop_qa = json.loads(row['crop_qa'])['qa_pairs'] 
            
            if (isinstance(crop_qa,list)):
                for c in crop_qa:
                    conversation = get_conversation(c['question'], c['answer'], crop_img, model_type)
                    json_data.append(conversation)
            else:
                conversation = get_conversation(crop_qa['question'], crop_qa['answer'], crop_img, model_type)
                json_data.append(conversation)
        return json_data


def save_train_jsonl(data, output_path, data_type, model_type, img_type,stage_num=1):
    output_dir_path = os.path.join(output_path, model_type)
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)
    if img_type not in ["transfer","mix"]:
        out = os.path.join(output_dir_path, f"stage{stage_num}_{data_type}_{img_type}.jsonl")
    else:
        out = os.path.join(output_dir_path, f"{stage_num}stage_{data_type}_{img_type}.jsonl")
    print("path",out)
    with open(out, 'w', encoding='utf-8') as f:
        for entry in data:
            json_line = json.dumps(entry)
            f.write(json_line + '\n')

def process_and_save(dataset, data_type, model_type):
    for data_path in dataset:
        if "bottom" not in data_path:
            top_file_data = data_path + '/top_raw_augment_qa_description_map_df.csv'
            side_anomaly_data = data_path + '/side_anomaly_augment_qa_description_map_df.csv'
            side_normally_data = data_path + '/side_normally_augment_qa_description_map_df.csv'

            top_stage1_data = construct_data(top_file_data, 'top', model_type)
            img_type="top"
            save_train_jsonl(top_stage1_data, data_path, data_type, model_type,img_type)

            side_anomaly_stage1_data = construct_data(side_anomaly_data, 'side_anomaly', model_type)
            img_type="side_anomaly"
            save_train_jsonl(side_anomaly_stage1_data, data_path, data_type, model_type,img_type)

            side_normally_stage1_data = construct_data(side_normally_data, 'side_normally', model_type)
            img_type="side_normally"
            save_train_jsonl(side_normally_stage1_data, data_path, data_type, model_type,img_type)
        else: 
            
            side_normally_data = data_path + '/mix_raw_augment_qa_description_map_df.csv'
            side_normally_stage1_data = construct_data(side_normally_data, 'mix', model_type)
            img_type="mix"
            save_train_jsonl(side_normally_stage1_data, data_path, data_type, model_type,img_type)
# split train test val
def merge_train_jsonl_files(root_dir, output_dir, data_type):
    pattern = os.path.join(root_dir, '**', f'*{"stage1"}*{data_type}*.jsonl')
    files_to_merge = glob.glob(pattern, recursive=True)

    if not files_to_merge:
        print("No matching JSONL files found.")
        return
    
    unique_lines_bottom = set()
    unique_lines_side = set()
    unique_lines_other = set()
    unique_lines_all = set()

    for file_path in files_to_merge:
        try:
            with open(file_path, 'r', encoding='utf-8') as f:
                for line in f:
                    stripped_line = line.strip()
                    if stripped_line:
                        try:
                            data = json.loads(stripped_line)
                            unique_lines_all.add(stripped_line + '\n')  # 添加到总集合

                            # 根据 images 中的内容分类
                            images = data.get('images', [])
                            if ("bottom" in images[0]):
                                unique_lines_bottom.add(stripped_line + '\n')
                            elif ("side"in images[0]):
                                unique_lines_side.add(stripped_line + '\n')
                            else:
                                unique_lines_other.add(stripped_line + '\n')
                        except json.JSONDecodeError:
                            print(f"Invalid JSON object in {file_path}: {stripped_line}")
        except OSError as e:
            print(f"Error reading {file_path}: {e}")

    # dir name
    output_files = {
        'all': output_dir+f"qwen_stage1_{data_type}.jsonl",
        'bottom': output_dir+f"bottom/qwen_stage1_{data_type}.jsonl",
        'side': output_dir+f"side/qwen_stage1_{data_type}.jsonl",
        'other':output_dir+f"top/qwen_stage1_{data_type}.jsonl",
    }

    # 将唯一行写入对应的输出文件
    for key, lines in [('all', unique_lines_all), ('bottom', unique_lines_bottom),
                       ('side', unique_lines_side), ('other', unique_lines_other)]:
        try:
            output_dir = os.path.dirname(output_files[key])
            os.makedirs(output_dir, exist_ok=True)  # 创建目录（如果不存在）
            with open(output_files[key], 'w', encoding='utf-8') as out_f:
                out_f.writelines(lines)
            print(f"Merged and categorized into {output_files[key]}")
        except OSError as e:
            print(f"Error writing to {output_files[key]}: {e}")

def to3(directory):
    """"create bottom view data"""
    # merging jsonl
    all_data = []
    for filename in os.listdir(directory):
        if filename.endswith('.jsonl') and "1stage" in filename:

            with open(os.path.join(directory, filename), 'r', encoding='utf-8') as file:
                for line in file:
                    all_data.append(json.loads(line))

    # delete old jsonl
    for filename in os.listdir(directory):
        if filename.endswith('.jsonl') and "1stage" in filename:
            os.remove(os.path.join(directory, filename))

    # shuffle data
    random.seed(42)  
    random.shuffle(all_data)


    total_items = len(all_data)
    train_split = int(0.7 * total_items)
    valid_split = train_split + int(0.2 * total_items)

  
    train_data = all_data[:train_split]
    valid_data = all_data[train_split:valid_split]
    test_data = all_data[valid_split:]

    def save_to_jsonl(data, filename):
        with open(filename, 'w', encoding='utf-8') as f:
            for item in data:
                f.write(json.dumps(item, ensure_ascii=False) + '\n')

    save_to_jsonl(train_data, os.path.join(directory, 'stage1_train.jsonl'))
    save_to_jsonl(valid_data, os.path.join(directory, 'stage1_valid.jsonl'))
    save_to_jsonl(test_data, os.path.join(directory, 'stage1_test.jsonl'))

if __name__=="__main__":
    
    data = "./Alphatool_data"
    from datetime import datetime
    now = datetime.now()
    formatted_time = now.strftime("%m%d-%H%M")
     
    data_list = [d for d in glob.glob(data + '/*') if os.path.isdir(d)   ]
    model_type = "qwen2vl"
    data_type=["train","val","test"]
    random.seed(42)  
    random.shuffle(data_list)  

    count=0
    total_length = len(data_list)
    train_ratio = 0.7
    val_ratio = 0.2
    test_ratio = 0.1
    count1=0
    train_split = int(total_length * train_ratio)
    val_split = train_split + int(total_length * val_ratio)
    
    train_data = data_list[:train_split]
    val_data = data_list[train_split:val_split]
    test_data = data_list[val_split:]
    print(len(train_data),train_data)
    print(len(val_data),val_data)
    print(len(test_data),test_data) 

    process_and_save(val_data, 'val', model_type)
    process_and_save(train_data, 'train', model_type)
    process_and_save(test_data, 'test', model_type)
    to3(os.path.join("./Alphatool_data/bottom",model_type))

    print(count)
    print(count1)
    root_directory = './Alphatool_data'

    data_types=["train","val","test"]
    for data_type in data_types:

        output_dir = f'./stage1_{formatted_time}/' 
        os.makedirs(output_dir, exist_ok=True)

        merge_train_jsonl_files(root_directory, output_dir,data_type)

