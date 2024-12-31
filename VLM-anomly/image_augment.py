import cv2
import numpy as np
import os
from glob import glob
import random
import pandas as pd
import shutil
from tqdm import tqdm
from PIL import Image

def augment_(img_path, output_folder, img_name):
    """
    Applies a series of random augmentations to an input image and returns the augmented image and its path.

    Parameters:
    img_path (str): The file path of the input image.
    output_folder (str): The folder where the augmented image will be saved.
    img_name (str): The name of the augmented image file.

    Returns:
    tuple: A tuple containing:
        - augmented_img (numpy.ndarray): The augmented image.
        - enhanced_img_path (str): The file path where the augmented image is saved.

    Functionality:
    1. Reads the input image using the provided `img_path`.
    2. Applies a random selection of augmentation techniques, including:
       - Rotate: Random rotation between -10 and 10 degrees.
       - Scale: Random scaling between 0.9 and 1.02 times.
       - Brightness: Adjust brightness by a random value in the range [-10, 10].
       - Color Jitter: Slight adjustments to saturation and brightness in the HSV color space.
       - Gamma: Adjusts gamma levels for brightness correction.
       - Contrast: Slightly modifies contrast.
       - Gaussian Noise: Adds Gaussian noise with a standard deviation of 0.3.
       - Salt and Pepper Noise: Introduces random salt-and-pepper noise.
       - Motion Blur: Simulates motion blur using a kernel of size 3.
    3. The selected augmentations are applied sequentially to the image.
    4. The augmented image is saved with a name reflecting the applied augmentations.

    Notes:
    - The function ensures the image is read properly. If the image cannot be read, the program exits with an error message.
    - The "Gaussian Noise" augmentation is always applied at the end in addition to the randomly chosen augmentations.

    Example:
    Suppose the input image is "input.jpg" located in the "images" folder. After running:
        augment_("images/input.jpg", "output", "augmented.jpg")
    The function will save the augmented image in the "output" folder with a filename indicating the applied augmentations.
    """
    img_path = os.path.abspath(img_path)
    img_pil = Image.open(img_path)  
    img = cv2.cvtColor(np.array(img_pil), cv2.COLOR_RGB2BGR) 

    if img is None:
        print(f"Unable to read image: {img_path}")
        exit()
        #continue

    # Randomly select a combination of augmentation operations
    augmentations = [
        "rotate", "scale", "brightness", "color_jitter", "gamma", "contrast",
        "gaussian_noise", "salt_pepper_noise", "motion_blur"
    ]
    num_augmentations = random.randint(1, 2)#len(augmentations))
    chosen_augmentations = random.sample(augmentations, num_augmentations)

    augmented_img = img.copy()
    height, width = img.shape[:2]

    for aug in chosen_augmentations+['gaussian_noise']:
        if aug == "rotate":
            # Random rotation between -10 and 10 degrees
            angle = random.uniform(-10, 10)
            center = (width // 2, height // 2)
            M_rot = cv2.getRotationMatrix2D(center, angle, 1.0)
            #augmented_img = cv2.warpAffine(augmented_img, M_rot, (width, height))
            augmented_img = cv2.warpAffine(augmented_img, M_rot, (width, height), flags=cv2.INTER_LINEAR, borderMode=cv2.BORDER_REFLECT_101)

        elif aug == "scale":
            # Random scaling (0.9 to 1.1 times)
            scale = random.uniform(0.9, 1.02)
            augmented_img = cv2.resize(augmented_img, None, fx=scale, fy=scale, interpolation=cv2.INTER_LINEAR)

        elif aug == "brightness":
            # Adjust brightness (random increase or decrease by up to 10 units)
            brightness_offset = random.randint(-10, 10)
            augmented_img = cv2.convertScaleAbs(augmented_img, alpha=1, beta=brightness_offset)

        elif aug == "color_jitter":
            # Apply slight color jitter
            hsv_img = cv2.cvtColor(augmented_img, cv2.COLOR_BGR2HSV)
            hsv_img = np.array(hsv_img, dtype=np.float64)
            hsv_img[..., 1] = hsv_img[..., 1] * random.uniform(0.9, 1.1)  # Saturation adjustment
            hsv_img[..., 2] = hsv_img[..., 2] * random.uniform(0.9, 1.1)  # Brightness adjustment
            hsv_img[..., 1:3] = np.clip(hsv_img[..., 1:3], 0, 255)
            augmented_img = cv2.cvtColor(np.array(hsv_img, dtype=np.uint8), cv2.COLOR_HSV2BGR)

        elif aug == "gamma":
            # Apply gamma correction (slight adjustment)
            gamma = random.uniform(0.9, 1.1)
            inv_gamma = 1.0 / gamma
            table = np.array([(i / 255.0) ** inv_gamma * 255 for i in np.arange(0, 256)]).astype(np.uint8)
            augmented_img = cv2.LUT(augmented_img, table)

        elif aug == "contrast":
            # Adjust contrast (slight adjustment)
            alpha = random.uniform(0.9, 1.1)  # Contrast factor
            augmented_img = cv2.convertScaleAbs(augmented_img, alpha=alpha, beta=0)

        elif aug == "gaussian_noise":
            # Add Gaussian noise
            gauss = np.random.normal(0, 0.3, augmented_img.shape).astype(np.uint8)
            augmented_img = cv2.add(augmented_img, gauss)

        elif aug == "salt_pepper_noise":
            # Add salt and pepper noise
            salt_pepper_ratio = 0.01
            num_salt = np.ceil(salt_pepper_ratio * augmented_img.size * 0.5)
            num_pepper = np.ceil(salt_pepper_ratio * augmented_img.size * 0.5)
            salt_pepper_img = augmented_img.copy()
            # Add salt noise
            coords = [np.random.randint(0, i - 1, int(num_salt)) for i in augmented_img.shape]
            salt_pepper_img[coords[0], coords[1], :] = 255
            # Add pepper noise
            coords = [np.random.randint(0, i - 1, int(num_pepper)) for i in augmented_img.shape]
            salt_pepper_img[coords[0], coords[1], :] = 0
            augmented_img = salt_pepper_img

        elif aug == "motion_blur":
            # Apply motion blur
            size = 3
            kernel_motion_blur = np.zeros((size, size))
            kernel_motion_blur[int((size - 1) / 2), :] = np.ones(size)
            kernel_motion_blur = kernel_motion_blur / size
            augmented_img = cv2.filter2D(augmented_img, -1, kernel_motion_blur)

    # Save the augmented image
    aug_str = '_'.join(chosen_augmentations)
    enhanced_img_path = os.path.join(output_folder, f"{aug_str}_{img_name}")
    
    return augmented_img, enhanced_img_path

def augment_images(raw_map_df, new_map_df):
    
    for _, raw_crop_img in tqdm(raw_map_df.iterrows(), total=len(raw_map_df)):
        raw_img_path = raw_crop_img['raw_img_path']
        crop_img_path = raw_crop_img['crop_img_path']
        img_class = raw_crop_img['img_class']

        raw_img_name = raw_img_path.split('/')[-1]
        crop_img_name = crop_img_path.split('/')[-1]
        output_dir = '/'.join(raw_img_path.split('/')[:-3])+'/augment'
        
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        augmented_img, enhanced_img_path = augment_(raw_img_path, output_dir, raw_img_name)
        # Pillow save
        if os.path.exists(enhanced_img_path):
            pass
        else:
            img_pil = Image.fromarray(cv2.cvtColor(augmented_img, cv2.COLOR_BGR2RGB))
            img_pil.save(enhanced_img_path, quality=100)

        augmented_crop_img, enhanced_crop_img_path = augment_(crop_img_path, output_dir, crop_img_name)
        # Pillow save
        if os.path.exists(enhanced_crop_img_path):
            pass
        else:
            img_pil = Image.fromarray(cv2.cvtColor(augmented_crop_img, cv2.COLOR_BGR2RGB))
            img_pil.save(enhanced_crop_img_path, quality=100)

        new_row = pd.DataFrame({"raw_img_path": [raw_img_path],'crop_image_path':[crop_img_path], 
                                "augment_raw_img_path":[enhanced_img_path],"augment_crop_img_path":[enhanced_crop_img_path],
                                'img_class': [img_class]})
        new_map_df = pd.concat([new_map_df, new_row], ignore_index=True)

    return new_map_df

