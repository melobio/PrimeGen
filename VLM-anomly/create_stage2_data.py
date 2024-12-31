import os
import random
import json
import pandas as pd
import glob

def get_conversation(query, response, img_path, model_type):
    data_type_list_1 = ['glm4v']
    data_type_list_2 = ['phi35', 'MiniCPM','qwen2vl','internvl2', 'llama32']
    data_type_list_3 = ['phi3']
    
    temp=img_path.replace("./Alphatool_data", "/home/Tli/workspace/PCR/Alphatool_data")

    #/home/YJhou/compute/VLM-2024
    if model_type in data_type_list_1:
        conversation= {"query": query, "response": response, "images": [temp]}
    elif model_type in data_type_list_2:
        conversation = {"query": f"<image>{query}", "response": response, "images": [temp]}
    elif model_type in data_type_list_3:
        conversation = {"conversations":[{"from": "user", "value": f"<img>{temp}</img>{query}"},{"from": "assistant", "value": response}]}
    else:
        print(f"error {model_type} not supported")

    return conversation


def construct_data(data_path,data_type, model_type):
    json_data=[]
    global count,count1
    df = pd.read_csv(data_path)
    if data_type == "top":
        for idx,row in df.iterrows():  
            temp=row['layout_qa'].replace("\'","\"")
             
            layout_qa = json.loads(temp, strict=False)  
            conversation = get_conversation(layout_qa['Q'][:-1]+"and return the description in JSON format just like {'slot 1': 'xxx', 'slot 2': 'xxx'}.", layout_qa['A'], row['augment_raw_img_path'], model_type)
            json_data.append(conversation)

            conversation = get_conversation(layout_qa['Q'][:-1]+"and return the description in JSON format just like {'slot 1': 'xxx', 'slot 2': 'xxx'}.", layout_qa['A'], row['raw_img_path'], model_type)
            json_data.append(conversation)

            temp=row['thermocycler_qa'].replace("\'","\"")
            a=json.loads(temp, strict=False) ["A"]
            thermocycler_q = r"Please determine the status of the AlphaTool thermal cycler during its operation and return the information in JSON format, where the value can only be none, open, or close, with the key being thermocycler_status."
            conversation = get_conversation(thermocycler_q,  a , row['augment_raw_img_path'], model_type)
            json_data.append(conversation)
                
        return json_data
    elif data_type == "side_anomaly":
        for idx,row in df.iterrows():  
            temp=row['stage2_qa'].replace("\'", "\"")
            qa_list = json.loads(temp) 
            conversation = get_conversation(qa_list['Q'][:-1]+"and return the description in JSON format.", qa_list['A'], row['augment_crop_img_path'], model_type)
            json_data.append(conversation) 
        return json_data

    elif data_type == "side_normally":
        for _,row in df.iterrows():  
            temp=row['stage2_qa'].replace("\'", "\"")
            qa_list = json.loads(temp) 
            conversation = get_conversation(qa_list['Q'][:-1]+"and return the description in JSON format.", qa_list['A'], row['augment_crop_img_path'], model_type)
            json_data.append(conversation) 
        return json_data
    else:
        for _,row in df.iterrows():  
            temp=row['stage2_qa'].replace("\'", "\"")
            qa_list = json.loads(temp) 
            conversation = get_conversation(qa_list['Q'][:-1]+"and return the description in JSON format.", qa_list['A'], row['augment_crop_img_path'], model_type)
            json_data.append(conversation) 
            #print(row)
            #exit()
            conversation = get_conversation(qa_list['Q'][:-1]+"and return the description in JSON format.", qa_list['A'], row['crop_image_path'], model_type)
            json_data.append(conversation)
        return json_data


def save_train_jsonl(data, output_path, data_type, model_type, img_type,stage_num=2):
    output_dir_path = os.path.join(output_path, model_type)
    if not os.path.exists(output_dir_path):
        os.makedirs(output_dir_path)
    if img_type not in ["transfer","mix"]:
        out = os.path.join(output_dir_path, f"stage{stage_num}_{data_type}_{img_type}.jsonl")
    else:
        out = os.path.join(output_dir_path, f"{stage_num}stage_{data_type}_{img_type}.jsonl")

    with open(out, 'w', encoding='utf-8') as f: 
        for entry in data:
            json_line = json.dumps(entry)
            f.write(json_line + '\n')

def process_and_save(dataset, data_type, model_type):
    for data_path in dataset:
        if "bottom" not in data_path:
            top_file_data = data_path + '/stage2_top_raw_augment_qa_description_map_df.csv'
            side_anomaly_data = data_path + '/stage2_side_anomaly_augment_qa_description_map_df.csv'
            
            side_normally_data = data_path + '/stage2_side_normally_augment_qa_description_map_df.csv'

            top_stage2_data = construct_data(top_file_data, 'top', model_type)
            img_type="top"
            save_train_jsonl(top_stage2_data, data_path, data_type, model_type,img_type)

            side_anomaly_stage2_data = construct_data(side_anomaly_data, 'side_anomaly', model_type)
            img_type="side_anomaly"
            save_train_jsonl(side_anomaly_stage2_data, data_path, data_type, model_type,img_type)

            side_normally_stage2_data = construct_data(side_normally_data, 'side_normally', model_type)
            img_type="side_normally"
            save_train_jsonl(side_normally_stage2_data, data_path, data_type, model_type,img_type)
        else: 
            # side_anomaly_data = data_path + '/stage2_transfer_augment_qa_description_map_df.csv'
            # side_anomaly_stage2_data = construct_data(side_anomaly_data, 'transfer', model_type)
            # img_type="transfer"
            # save_train_jsonl(side_anomaly_stage2_data, data_path, data_type, model_type,img_type)
            
            side_normally_data = data_path + '/stage2_mix_augment_qa_description_map_df.csv'
            side_normally_stage2_data = construct_data(side_normally_data, 'mix', model_type)
            img_type="mix"
            save_train_jsonl(side_normally_stage2_data, data_path, data_type, model_type,img_type)
            
def merge_train_jsonl_files(root_dir, output_dir, data_type):
    # 使用 glob 查找所有包含 'stage2' 和 data_type 的 .jsonl 文件
    pattern = os.path.join(root_dir, '**', f'*{"stage2"}*{data_type}*.jsonl')
    files_to_merge = glob.glob(pattern, recursive=True)

    if not files_to_merge:
        print("No matching JSONL files found.")
        return
    
    # 初始化集合用于去重和分类存储
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

    # 定义输出文件名
    output_files = {
        'all': output_dir+f"qwen_stage2_{data_type}.jsonl",
        'bottom': output_dir+f"bottom/qwen_stage2_{data_type}.jsonl",
        'side': output_dir+f"side/qwen_stage2_{data_type}.jsonl",
        'other':output_dir+f"top/qwen_stage2_{data_type}.jsonl",
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
    # 读取并合并所有jsonl文件
    all_data = []
    for filename in os.listdir(directory):
        if filename.endswith('.jsonl') and "2stage" in filename:
            with open(os.path.join(directory, filename), 'r', encoding='utf-8') as file:
                for line in file:
                    all_data.append(json.loads(line))

    # 删除原来的jsonl文件
    # for filename in os.listdir(directory): 
    #     if filename.endswith('.jsonl') and "2stage" in filename:
    #         print("remove",filename)
    #         os.remove(os.path.join(directory, filename))

    # 随机打乱数据以确保随机性
    random.seed(42)  # 确保结果的可重复性
    random.shuffle(all_data)

    # 计算分割点
    total_items = len(all_data)
    train_split = int(0.7 * total_items)
    valid_split = train_split + int(0.2 * total_items)

    # 分割数据集
    train_data = all_data[:train_split]
    valid_data = all_data[train_split:valid_split]
    test_data = all_data[valid_split:]

    # 将数据保存到对应的json文件中
    def save_to_jsonl(data, filename):
        with open(filename, 'w', encoding='utf-8') as f:
            for item in data:
                f.write(json.dumps(item, ensure_ascii=False) + '\n')

    save_to_jsonl(train_data, os.path.join(directory, 'stage2_train.jsonl'))
    save_to_jsonl(valid_data, os.path.join(directory, 'stage2_valid.jsonl'))
    
    save_to_jsonl(test_data, os.path.join(directory, 'stage2_test.jsonl'))
def check_image_files(jsonl_file_path, base_dir=''):
    """
    Read the jsonl file and check if the image files corresponding to each img_path exist.
    
    :param jsonl_file_path: The path to the jsonl file containing the img_path key.
    :param base_dir: The base directory for the image files, default is an empty string, indicating a relative path.
    """
    # 记录不存在的文件
    missing_files = []

    with open(jsonl_file_path, 'r', encoding='utf-8') as file:
        for line in file:
            try:
                data = json.loads(line.strip())
                img_path = data.get('images', None)[0]

                if img_path is None:
                    print("Warning: img_path key not found in one of the lines.")
                    continue


                full_img_path =  img_path

  
                if not os.path.exists(full_img_path):
                    missing_files.append(full_img_path) 

            except json.JSONDecodeError as e:
                print(f"Failed to decode JSON line: {e}")
                continue

    if missing_files:
        print("\nSummary of missing files:")
        for mf in missing_files:
            print(mf)
    else:
        print("All image files exist.")

if __name__=="__main__": 
    data = "./Alphatool_data"
    data_list = [d for d in glob.glob(data + '/*') if os.path.isdir(d)  ]

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
    to3("./Alphatool_data/bottom/qwen2vl")
    from datetime import datetime


    now = datetime.now()
    formatted_time = now.strftime("%m%d-%H%M")

    print(count)
    root_directory = './Alphatool_data'
    data_types=["train","val","test"]
    for data_type in data_types:
        output_dir = f'./stage2_{formatted_time}/' 
        os.makedirs(output_dir, exist_ok=True)

        merge_train_jsonl_files(root_directory, output_dir,data_type)
