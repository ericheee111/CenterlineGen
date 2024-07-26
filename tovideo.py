import cv2
import os

def make_video(image_folder, video_name, fps):
    images = [img for img in os.listdir(image_folder) if img.endswith(".png")]
    # Ensure the images are in the correct order
    images.sort(key=lambda x: int(x.split('out')[1].split('.')[0]))

    # Extract the first image to get the dimensions
    frame = cv2.imread(os.path.join(image_folder, images[0]))
    height, width, layers = frame.shape

    # Define the codec and create VideoWriter object
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # 'mp4v' is the codec for .mp4 format
    video = cv2.VideoWriter(video_name, fourcc, fps, (width, height))
    
    # Add images to video
    for i, image in enumerate(images):
        img = cv2.imread(os.path.join(image_folder, image))
        cv2.putText(img, str(i+1), (10, 30), cv2.FONT_HERSHEY_SIMPLEX, 1, (0, 0, 255), 2)
        video.write(img)

    # Release the video
    video.release()
    print("Video creation complete.")

# Usage
image_folder = "C:/Users/Eric/Desktop/study/gacrnd/CenterlineGen/CenterlineGen"
video_name = 'output_video.mp4'       # Name of the output video file
fps = 10                              # Frames per second

make_video(image_folder, video_name, fps)