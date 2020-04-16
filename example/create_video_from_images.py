import imageio
import PIL.Image
import PIL.ImageOps
import PIL.ImageDraw
import os
from path import Path
from math import floor


# -- To save images from PlantGL
# from openalea.plantgl.all import Scene as pglScene, Viewer
# c_scene_sky.plot() # =PARa_sky, minval=0., maxval=700.)
# Viewer.camera.setPosition([0, 1, 0])
# Viewer.frameGL.saveImage('video/scene_{}.png'.format(t))


# -- Definition of a function that can resize a list of images and make a movie from it:
def cropping_resizing_and_film_making(cropping=False, resizing=False, framing=True, cropping_area=None, resizing_size=None, film_making=True,
                                      original_images_directory='video', film_directory='video', sampling_frequency=1):
    """
    This function can resize a list of images contained in the directory "video", store them in a new directory, and create a .gif movie from it.
    :param bool cropping: a Boolean allowing to crop
    :param bool resizing: a Boolean allowing to resize
    :param bool framing: a Boolean allowing to add a frame in which the simulation time is displayed
    :param tuple cropping_area: Cropping area  (x1, y1, x2, y2)
    :param tuple resizing_size: a tupple containing the width and the height of the resized image
    :param bool film_making: a Boolean allowing to create a video file "root_movie.gif"
    :param str original_images_directory: the name of the directory where the original images will be read to be cropped or resized
    :param str film_directory: the name of the directory where the images (cropped, resized or original) will be read to create the movie
    :param int sampling_frequency: the code will sample the images at the sampling_frequency (to reduce the size of the movie and/or accelerate the movie)
    :return:
    """

    # Getting a list of the names of the images found in the directory "video":
    filenames = Path(original_images_directory).glob('*.png')
    filenames = sorted(filenames)

    # We define the final number of images that will be considered, based on the "sampling_frequency" variable:
    number_of_images = floor(len(filenames) / float(sampling_frequency)) + 1

    # We create a new directory "video_adjusted" that will contain the cropped and/or resized images, if necessary:
    if resizing or cropping:
        video_dir = 'video_adjusted'
        # If this directory doesn't exist:
        if not os.path.exists(video_dir):
            # Then we create it:
            os.mkdir(video_dir)
        else:
            # Otherwise, we delete all the images that are already present inside:
            for root, dirs, files in os.walk(video_dir):
                for file in files:
                    os.remove(os.path.join(root, file))

    # 1. CROPPING THE IMAGES:
    if cropping:
        # We resize the images:
        print "Cropping the images..."
        number = 0
        count = 0
        remaining_images = number_of_images

        # We cover each image in the directory:
        for filename in filenames:
            # The count is increased:
            count += 1
            # If it corresponds to the target number, the image is added to the gif:
            if count == sampling_frequency:
                print "Please wait:", str(int(remaining_images)), "image(s) left"

                simulation_time = int(filename[:-4][-4:])

                im = PIL.Image.open(filename)
                im_cropped = im.crop(cropping_area)
                image_name = os.path.join(video_dir, 'adjusted_image_%.4d.png')
                im_cropped.save(image_name % simulation_time)
                remaining_images = remaining_images - 1
                # We reset the count to 0:
                count = 0
        print "The images have been resized!"

    if framing:
        if cropping:
            # Getting a list of the names of the images found in the directory "video":
            filenames = Path('video_adjusted').glob('*.png')
            filenames = sorted(filenames)

        # We add a frame to the images:
        print "Adding frame to the images..."
        count = 0
        remaining_images = number_of_images

        # We cover each image in the directory:
        for filename in filenames:
            # The count is increased:
            count += 1
            # If it corresponds to the target number, the image is added to the gif:
            if count == sampling_frequency:
                print "Please wait:", str(int(remaining_images)), "image(s) left"

                simulation_time = int(filename[:-4][-4:])

                im = PIL.Image.open(filename)
                # Frame
                im_framed = PIL.ImageOps.expand(im, border=(0, 0, 0, 50), fill='white')
                # Display simulation time
                width, height = im.size
                im_framed2 = PIL.ImageDraw.Draw(im)
                im_framed2.text((10, height - 10), 't = %.4d h' % simulation_time, fill=(255,255,255))
                # Save
                image_name = os.path.join(video_dir, 'adjusted_image_%.4d.png')
                im_framed.save(image_name % simulation_time)
                remaining_images = remaining_images - 1
                # We reset the count to 0:
                count = 0
        print "The frames have been added!"

    # 3. COMPRESSING THE IMAGES:
    if resizing:
        if framing or cropping:
            # Getting a list of the names of the images found in the directory "video":
            filenames = Path('video_adjusted').glob('*.png')
            filenames = sorted(filenames)
        # We resize the images:
        print "Resizing the images..."
        number = 0
        count = 0
        remaining_images = number_of_images

        # We cover each image in the directory:
        for filename in filenames:
            # The count is increased:
            count += 1
            # If it corresponds to the target number, the image is added to the gif:
            if count == sampling_frequency:
                print "Please wait:", str(int(remaining_images)), "image(s) left"
                im = PIL.Image.open(filename)
                im_resized = im.resize(resizing_size, resample=0)
                image_name = os.path.join(video_dir, 'adjusted_image_%.4d.png')
                im_resized.save(image_name % number, quality=20, optimize=True)
                number = number + 1
                remaining_images = remaining_images - 1
                # We reset the count to 0:
                count = 0
        print "The images have been resized!"

    # 4. CREATING THE VIDEO FILE:
    if film_making:

        print "Making the video..."

        filenames = Path(film_directory).glob('*.png')
        filenames = sorted(filenames)

        with imageio.get_writer('movie.gif', mode='I', fps=24) as writer:
            if resizing:
                sampling_frequency = 1
            else:
                sampling_frequency = sampling_frequency
            remaining_images = floor(len(filenames) / float(sampling_frequency)) + 1
            print remaining_images, "images are considered at this stage."
            # We add the first image:
            filename = filenames[0]
            image = imageio.imread(str(filename))
            writer.append_data(image)
            # We reduce the number of images left:
            remaining_images = remaining_images - 1
            # We start the count at 0:
            count = 0
            # We cover each image in the directory:
            for filename in filenames:
                # The count is increased:
                count += 1
                # If it corresponds to the target number, the image is added to the gif:
                if count == sampling_frequency:
                    print "Please wait:", str(int(remaining_images)), "image(s) left"
                    image = imageio.imread(str(filename))
                    writer.append_data(image)
                    remaining_images = remaining_images - 1
                    # We reset the count to 0:
                    count = 0
        print "The video has been made!"


if __name__ == '__main__':
    cropping_resizing_and_film_making(cropping=True, framing=True, cropping_area=(200, 110, 550, 440), film_making=False,
                                      original_images_directory='video', film_directory='video_adjusted', sampling_frequency=500)
