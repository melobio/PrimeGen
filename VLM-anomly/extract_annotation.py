import xml.etree.ElementTree as ET
import os
import json

#Extract data based on annotations made using CVAT, for subsequent data construction.
def get_annotation(path, img_dir):
    img_label_list = []
    tree = ET.parse(path)
    root = tree.getroot()
    # Traverse all image nodes
    for image in root.findall('image'):
        image_id = image.get('id')
        image_name = image.get('name')
        width = image.get('width')
        height = image.get('height')
        print(f"Image ID: {image_id}, Name: {image_name}, Width: {width}, Height: {height}")
        
        labels = []
        coordinates = []
        
        # Iterate through all annotation boxes in the image
        for box in image.findall('box'):
            label = box.get('label')
            xtl = box.get('xtl')
            ytl = box.get('ytl')
            xbr = box.get('xbr')
            ybr = box.get('ybr')
            print(f"  Label: {label}, Coordinates: ({xtl}, {ytl}), ({xbr}, {ybr})")
            labels.append(label)
            coordinates.append([xtl, ytl, xbr, ybr])
        if labels:
            img_label_list.append({
                "img_path": f"{img_dir}/{image_name}",
                "labels": labels,
                "coordinates": coordinates
            })
        
    return img_label_list

def save_to_jsonl(data, jsonl_file):
    with open(jsonl_file, 'a', encoding='utf-8') as f:
        for entry in data:
            json.dump(entry, f, ensure_ascii=False)
            f.write('\n')


def check_if_processed(annotation_path, processed_files_log):
    # Check if the file has already been processed
    if not os.path.exists(processed_files_log):
        return False
    with open(processed_files_log, 'r', encoding='utf-8') as f:
        processed_files = f.read().splitlines()
    return annotation_path in processed_files

def log_processed_file(annotation_path, processed_files_log):
    # Log the path of the processed file
    with open(processed_files_log, 'a', encoding='utf-8') as f:
        f.write(f"{annotation_path}\n")

def main(annotation_path_list, jsonl_file, processed_files_log):
    img_label=[]
    for annotation_path in annotation_path_list:
        # Check if the file has already been handled
        if check_if_processed(annotation_path, processed_files_log):
            print(f"Skipping already processed file: {annotation_path}")
            continue
        print(annotation_path)
        data_name = '-'.join(annotation_path.split('/')[-2].split("-")[:-1])
        img_dir = '/'.join(annotation_path.split('/')[:-2]) +'/'+ data_name
        img_label_list=get_annotation(annotation_path, img_dir)
        img_label+=img_label_list
        # Save to a JSONL file
        save_to_jsonl(img_label_list, jsonl_file)
        # Record the processed files
        log_processed_file(annotation_path, processed_files_log)

if __name__=="__main__":
    annotation_path_list = ['OT2-INCREMENT-DATA-9-9-Annotation/annotations.xml','OT2-INCREMENT-DATA-Annotation/annotations.xml', 'OT2-DATA-Annotation/annotations.xml'] 
    jsonl_file = 'annotations_data.jsonl'
    processed_files_log = 'processed_files.log'
    main(annotation_path_list, jsonl_file, processed_files_log)

