# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 20:28:12 2024

@author: 20481756
"""

from PIL import Image
format_name = 'png'
image1 = Image.open('../IllustrationFigure/mechanism.' + format_name)
image2 = Image.open('../Figure/Figure 3de.' + format_name)
width1, height1 = image1.size
width2, height2 = image2.size
aspect_ratio = height1 / width1
new_height1 = int(width2 * aspect_ratio)
resized_image1 = image1.resize((width2, new_height1))
combined_height = new_height1 + height2
merged_image = Image.new('RGB', (width2, combined_height))
merged_image.paste(resized_image1, (0, 0))
merged_image.paste(image2, (0, new_height1))
merged_image.save('../Figure/Figure 3.' + format_name)