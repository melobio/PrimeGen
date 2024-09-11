import os 
import random
import json
import pandas as pd

#Here, we are constructing various types of VLM data, with the underlying format compatible with Swift training data.
def construct_data(data_dir,output_path):
    #TODO: Adaptively construct various types of data
    seed = random.seed(42)
    df = pd.read_csv(data_dir)

    #json_data_glm = []
    json_data_phi3 = []
    #json_data_qwen2vl = []
    for index, row in df.iterrows():    
        #phi3
        conversation_phi3={"conversations":[{"from": "user", "value": f"<img>/home/YJhou/compute/VLM-2024/PCR/Data/Data_process{row['img_path']}</img>{row['Q']}"},{"from": "assistant", "value": row["A"]}]}
        
        #glm
        #conversation_glm ={"query": row["Question"], "response": row["Answer"], "images": [row["Path"]]}
        
        #qwen2vl
        #conversation_qwen2vl={"query":f"<image>{row['Q']}","response":row["A"],"images":[f"/home/YJhou/compute/VLM-2024/PCR/Data/Data_process{row['img_path']}"]}
        

        #json_data_glm.append(conversation_glm)
        json_data_phi3.append(conversation_phi3)
        #json_data_qwen2vl.append(conversation_qwen2vl)

        # 随机打乱
        # random.shuffle(json_data_glm)
        random.shuffle(json_data_phi3)
        #random.shuffle(json_data_qwen2vl)

        #length = len(json_data_glm)
        length = len(json_data_phi3)
        #length = len(json_data_qwen2vl)
        

        output_file_path_qwen2vl = os.path.join(output_path,'Stage_2_data_qwen2vl')
        os.makedirs(output_file_path_qwen2vl,exist_ok=True)

        #output_file_path_phi3 = os.path.join(output_path,'Stage_2_data_phi3')
        #os.makedirs(output_file_path_phi3,exist_ok=True)

        #output_file_path_glm = os.path.join(output_path,'Stage_2_data_glm')
        #os.makedirs(output_file_path_glm,exist_ok=True)
        
        train_num = int(length * 0.8)
        val_num = int(length * 0.1)

        #json_data_glm_train = json_data_glm[:train_num]
        #json_data_glm_val = json_data_glm[train_num:val_num+train_num]
        #json_data_glm_test = json_data_glm[val_num+train_num:]

        #json_data_phi3_train = json_data_phi3[:train_num]
        #json_data_phi3_val = json_data_phi3[train_num:val_num+train_num]
        #json_data_phi3_test = json_data_phi3[val_num+train_num:]
        
        json_data_qwen2vl_train = json_data_qwen2vl[:train_num]
        json_data_qwen2vl_val = json_data_qwen2vl[train_num:val_num+train_num]
        json_data_qwen2vl_test = json_data_qwen2vl[val_num+train_num:]
        
        for split in ["train", "val", "test"]:
            
            #data_phi3 = {
            #    "train": json_data_phi3_train,
            #    "val": json_data_phi3_val,
            #    "test": json_data_phi3_test
            #}[split]
            
            #data_glm = {
            #    "train": json_data_glm_train,
            #    "val": json_data_glm_val,
            #    "test": json_data_glm_test
            #}[split]

            data_qwen2vl = {
                  "train": json_data_qwen2vl_train,
                  "val": json_data_qwen2vl_val,
                  "test": json_data_qwen2vl_test
                  }[split]
            
            out_qwen2vl = os.path.join(output_file_path_qwen2vl, f"{split}.jsonl")
            with open(out_qwen2vl, 'w', encoding='utf-8') as f:
                for entry in data_qwen2vl:
                    json_line = json.dumps(entry)
                    f.write(json_line + '\n')

            #out_phi = os.path.join(output_file_path_phi3, f"{split}.jsonl")
            #with open(out_phi, 'w', encoding='utf-8') as f:
            #    for entry in data_phi3:
            #        json_line = json.dumps(entry)
            #        f.write(json_line + '\n')
    
            #out_glm = os.path.join(output_file_path_glm, f"{split}.jsonl")
            #with open(out_glm, 'w', encoding='utf-8') as f:
            #    for entry in data_glm:
            #        json_line = json.dumps(entry)
            #        f.write(json_line + '\n')

if __name__=="__main__":
    out = "./"
    data = "../Data_augment/Global_VQA.csv"
    construct_data(data,out)
