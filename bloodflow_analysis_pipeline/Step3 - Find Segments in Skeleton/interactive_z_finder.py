#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Interactive 3D Point Finder
This program allows users to click on a maximum intensity projection image
to automatically find the nearest point with value 1 in 3D space
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Button
import tkinter as tk
from tkinter import filedialog
from PIL import Image

def find_nearest_nonzero_3d(img_stack, x, y, search_radius):
    """Find the nearest point with value 1 near the Z-axis line of the clicked point"""
    height, width, depth = img_stack.shape
    best_x = x
    best_y = y
    best_z = 0
    best_value = 0
    min_distance = float('inf')
    
    # For each Z level, find the nearest point with value 1 in XY plane
    for z in range(depth):
        local_min_distance = float('inf')
        local_best_x = x
        local_best_y = y
        
        # Search in XY plane around the clicked point
        for dy in range(-search_radius, search_radius + 1):
            for dx in range(-search_radius, search_radius + 1):
                test_x = x + dx
                test_y = y + dy
                
                # Check if within image bounds
                if (0 <= test_x < width and 
                    0 <= test_y < height):
                    
                    test_value = img_stack[test_y, test_x, z]
                    
                    # If found a point with value 1
                    if test_value == 1:
                        # Calculate 2D distance in current XY plane
                        distance_2d = np.sqrt(dx**2 + dy**2)
                        if distance_2d < local_min_distance:
                            local_min_distance = distance_2d
                            local_best_x = test_x
                            local_best_y = test_y
        
        # If found a point in this Z level
        if local_min_distance < float('inf'):
            # Calculate 3D distance including Z difference from middle of stack
            z_weight = 0.01  # Adjust this to change importance of Z distance
            distance_3d = np.sqrt(local_min_distance**2 + (z_weight * (z-depth/2)**2))
            
            # Update best point if this is the closest so far
            if distance_3d < min_distance:
                min_distance = distance_3d
                best_x = local_best_x
                best_y = local_best_y
                best_z = z
                best_value = 1
    
    return best_x, best_y, best_z, min_distance, best_value

def read_tiff_stack(filepath):
    """Read TIFF image stack"""
    img = Image.open(filepath)
    images = []
    try:
        while True:
            images.append(np.array(img))
            img.seek(img.tell() + 1)
    except EOFError:
        pass
    
    return np.array(images).transpose(1, 2, 0)  # Convert to (H, W, D) format

class InteractiveZFinder:
    def __init__(self):
        self.fig, self.ax = plt.subplots(figsize=(10, 8))
        self.fig.canvas.manager.set_window_title('3D Point Finder')
        
        # Select file
        root = tk.Tk()
        root.withdraw()  # Hide main window
        self.filepath = filedialog.askopenfilename(
            title='Select TIFF File',
            filetypes=[('TIFF files', '*.tif;*.tiff')]
        )
        
        if not self.filepath:
            plt.close()
            return
            
        # Read image data
        print(f'Reading file: {self.filepath}')
        self.img_stack = read_tiff_stack(self.filepath)
        print(f'Image size: {self.img_stack.shape}')
        
        # Calculate maximum intensity projection
        self.max_projection = np.max(self.img_stack, axis=2)
        
        # Display maximum intensity projection
        self.im = self.ax.imshow(self.max_projection, cmap='summer')
        self.ax.set_title('Click anywhere to find the nearest point with value 1')
        
        # Add click event handler
        self.fig.canvas.mpl_connect('button_press_event', self.on_click)
        
        # Add close button
        self.close_ax = plt.axes([0.8, 0.02, 0.1, 0.04])
        self.close_button = Button(self.close_ax, 'Close')
        self.close_button.on_clicked(self.close)
        
        # Store references to markers and texts
        self.markers = []
        self.texts = []
        
        plt.show()
    
    def clear_markers(self):
        """Clear previous markers"""
        for marker in self.markers:
            marker.remove()
        for text in self.texts:
            text.remove()
        self.markers = []
        self.texts = []
    
    def on_click(self, event):
        """Handle click events"""
        if event.inaxes != self.ax:
            return
            
        # Clear previous markers
        self.clear_markers()
        
        # Get click coordinates
        x, y = int(event.xdata), int(event.ydata)
        
        # Find nearest point with value 1 in 3D space
        search_radius = 20
        best_x, best_y, best_z, min_distance, best_value = find_nearest_nonzero_3d(
            self.img_stack, x, y, search_radius
        )
        
        # Verify the found point's value
        if best_value > 0:
            actual_value = self.img_stack[best_y, best_x, best_z]
            print(f'Verifying point ({best_x},{best_y},{best_z}) value: {actual_value}')
            
            if actual_value != 1:
                print('Warning: Found point value is not 1!')
            
            # Mark original click position
            marker = self.ax.plot(x, y, 'x', color='gray', markersize=12, linewidth=2)[0]
            self.markers.append(marker)
            
            # If found a different position
            if best_x != x or best_y != y:
                # Draw connection line
                line = self.ax.plot([x, best_x], [y, best_y], '--', 
                                  color='gray', linewidth=1.5)[0]
                self.markers.append(line)
            
            # Mark found point
            marker = self.ax.plot(best_x, best_y, '+', color='red', 
                                markersize=15, linewidth=2)[0]
            self.markers.append(marker)
            marker = self.ax.plot(best_x, best_y, 'o', color='red', 
                                markersize=8, linewidth=2)[0]
            self.markers.append(marker)
            
            # Display coordinate information
            text_color = 'red' if actual_value == 1 else 'orange'
            if best_x != x or best_y != y:
                text = self.ax.text(best_x+10, best_y-10, 
                    f'3D Coord: ({best_x},{best_y},{best_z})\n'
                    f'Value={actual_value}\n3D Dist={min_distance:.1f}',
                    color=text_color, fontsize=10, fontweight='bold',
                    bbox=dict(facecolor='white', alpha=1))
            else:
                text = self.ax.text(best_x+10, best_y-10,
                    f'3D Coord: ({best_x},{best_y},{best_z})\nValue={actual_value}',
                    color=text_color, fontsize=10, fontweight='bold',
                    bbox=dict(facecolor='white', alpha=1))
            self.texts.append(text)
            
            # Display complete information in command line
            print(f'Found 3D point: ({best_x},{best_y},{best_z}), '
                  f'Value={actual_value}, Distance={min_distance:.1f}')
        else:
            # If no point with value 1 found
            marker = self.ax.plot(x, y, 'rx', markersize=12, linewidth=2)[0]
            self.markers.append(marker)
            text = self.ax.text(x+10, y-10, 'No point with value 1 found',
                color='red', fontsize=10, fontweight='bold',
                bbox=dict(facecolor='white', alpha=1))
            self.texts.append(text)
        
        self.fig.canvas.draw_idle()
    
    def close(self, event):
        """Close window"""
        plt.close(self.fig)

if __name__ == '__main__':
    app = InteractiveZFinder() 