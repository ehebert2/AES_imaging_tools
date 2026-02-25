import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev
from skimage.measure import profile_line
import tifffile
import os

def load_tiff_stack(filename, start_frame=0, num_frames=None):
    """
    Load TIFF image stack with support for selecting time range
    
    Parameters:
    filename : str
        Path to TIFF file
    start_frame : int
        Starting frame index (0-based)
    num_frames : int or None
        Number of frames to read, None means read to the end
        
    Returns:
    numpy.ndarray
        Image stack with shape (frames, height, width)
    """
    with tifffile.TiffFile(filename) as tif:
        total_frames = len(tif.pages)
        print(f"TIFF file information:")
        print(f"Total frames: {total_frames}")
        print(f"Image size: {tif.pages[0].shape}")
        
        # Calculate the actual frame range to read
        if start_frame >= total_frames:
            raise ValueError(f"Starting frame {start_frame} exceeds total frames {total_frames}")
        
        if num_frames is None:
            end_frame = total_frames
        else:
            end_frame = min(start_frame + num_frames, total_frames)
        
        actual_frames = end_frame - start_frame
        print(f"Will read {actual_frames} frames, starting from frame {start_frame}")
        
        # Read frames in the selected range
        stack = []
        for i in range(start_frame, end_frame):
            stack.append(tif.pages[i].asarray())
        
        # Convert to numpy array
        stack = np.array(stack)
        print(f"Reading complete, image stack shape: {stack.shape}")
        return stack

def interpolate_curve(points, num_points=100):
    """
    Interpolate user-selected points to generate a smooth curve
    
    Parameters:
    points : numpy.ndarray
        Coordinates of user-selected points, shape (N, 2)
    num_points : int
        Number of points to generate after interpolation
        
    Returns:
    numpy.ndarray
        Interpolated point coordinates, shape (num_points, 2)
    """
    # Separate x and y coordinates
    x = points[:, 0]
    y = points[:, 1]
    
    # Use splprep for spline interpolation
    tck, u = splprep([x, y], s=0, k=min(3, len(points)-1))
    
    # Generate uniformly distributed parameter values
    u_new = np.linspace(0, 1, num_points)
    
    # Calculate interpolated points
    x_new, y_new = splev(u_new, tck)
    
    return np.column_stack((x_new, y_new))

def create_vessel_kymograph(image_stack, vessel_points, line_width=3, spacing=1):
    """
    Generate spatiotemporal map (kymograph) of the vessel
    
    Parameters:
    image_stack : numpy.ndarray
        3D image stack (time, height, width) or 2D image (height, width)
    vessel_points : numpy.ndarray
        Coordinate points of vessel centerline, shape (N, 2)
    line_width : int
        Width of vessel line (in pixels), must be odd
    spacing : float
        Spacing between interpolated points (in pixels)
        
    Returns:
    numpy.ndarray
        Generated kymograph with shape (time, profile_length)
    """
    if line_width % 2 == 0:
        raise ValueError("line_width must be odd")
        
    # Ensure image_stack is 3D
    if len(image_stack.shape) == 2:
        print("Warning: Input is a single frame image")
        image_stack = image_stack[np.newaxis, :, :]
    
    num_frames = image_stack.shape[0]
    print(f"Processing {num_frames} frames")
    print(f"Vessel line width: {line_width} pixels")
    
    # Calculate total length of the curve
    diff = np.diff(vessel_points, axis=0)
    segment_lengths = np.sqrt(np.sum(diff**2, axis=1))
    total_length = np.sum(segment_lengths)
    
    # Determine number of interpolation points based on total length
    num_points = int(total_length / spacing)
    print(f"Total curve length: {total_length:.1f} pixels")
    print(f"Number of sampling points: {num_points}")
    
    # Interpolate the curve
    interpolated_points = interpolate_curve(vessel_points, num_points)
    
    # Initialize kymograph matrix
    kymograph = np.zeros((num_frames, num_points-1))
    
    # Process each frame
    for t in range(num_frames):
        # Process each adjacent point pair on the curve
        for i in range(num_points-1):
            src = interpolated_points[i]
            dst = interpolated_points[i+1]
            
            # Swap x,y coordinates to match profile_line requirements
            src_yx = (src[1], src[0])
            dst_yx = (dst[1], dst[0])
            
            # Calculate profile for this segment
            profile = profile_line(image_stack[t], src_yx, dst_yx,
                                 mode='reflect',
                                 linewidth=line_width,
                                 reduce_func=np.mean)
            
            # Take the mean value of this segment as the intensity at this point
            kymograph[t, i] = np.mean(profile)
    
    return kymograph

def plot_kymograph(kymograph, start_frame, save_path=None):
    """
    Display and save kymograph
    
    Parameters:
    kymograph : numpy.ndarray
        Kymograph to display
    start_frame : int
        Starting frame index for correct time axis display
    save_path : str, optional
        Path to save the image
    """
    plt.figure(figsize=(12, 6))
    if len(kymograph.shape) == 1:
        kymograph = kymograph[np.newaxis, :]
    
    # 创建时间轴标签
    num_frames = kymograph.shape[0]
    plt.imshow(kymograph, aspect='auto', cmap='gray')
    plt.colorbar(label='Average Intensity')
    
    # 设置y轴（时间轴）的刻度
    frame_ticks = np.linspace(0, num_frames-1, 5, dtype=int)
    frame_labels = frame_ticks + start_frame
    plt.yticks(frame_ticks, frame_labels)
    
    plt.xlabel('Position along vessel (pixels)')
    plt.ylabel('Frame number')
    plt.title('Vessel Kymograph (Curved Path)')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def plot_vessel_path(image, points, interpolated_points=None, save_path=None):
    """
    Display vessel path on image with optional saving
    
    Parameters:
    image : numpy.ndarray
        Image to display
    points : numpy.ndarray
        User-selected points
    interpolated_points : numpy.ndarray, optional
        Interpolated points
    save_path : str, optional
        Save path (without extension)
    """
    plt.figure(figsize=(10, 10))
    plt.imshow(image, cmap='gray')
    
    # 绘制用户选择的点
    plt.plot(points[:, 0], points[:, 1], 'ro', label='Selected points')
    
    # 如果有插值点，也画出来
    if interpolated_points is not None:
        plt.plot(interpolated_points[:, 0], interpolated_points[:, 1], 'y-', 
                label='Interpolated path', linewidth=1)
    
    plt.legend()
    plt.title('Vessel Path')
    
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()

def save_vessel_path(points, interpolated_points, base_path):
    """
    Save vessel path data
    
    Parameters:
    points : numpy.ndarray
        User-selected points
    interpolated_points : numpy.ndarray
        Interpolated points
    base_path : str
        Base save path (without extension)
    """
    # Save original points and interpolated points
    path_data = {
        'original_points': points,
        'interpolated_points': interpolated_points
    }
    np.savez(f"{base_path}_path.npz", **path_data)
    print(f"- {base_path}_path.npz (vessel path data)")

def save_kymograph(kymograph, base_path):
    """
    Save kymograph in multiple formats
    
    Parameters:
    kymograph : numpy.ndarray
        Kymograph data to save
    base_path : str
        Base file path (without extension)
    """
    try:
        # Print information about input data
        print(f"\nSaving kymograph information:")
        print(f"Data shape: {kymograph.shape}")
        print(f"Data type: {kymograph.dtype}")
        print(f"Value range: [{kymograph.min():.2f}, {kymograph.max():.2f}]")
        
        # 1. Save as numpy array
        np_path = f"{base_path}.npy"
        np.save(np_path, kymograph)
        print(f"Successfully saved NPY file: {np_path}")
        
        # 2. Save as TIFF file
        tiff_path = f"{base_path}.tif"
        
        # Check if data contains NaN or inf
        if np.any(np.isnan(kymograph)) or np.any(np.isinf(kymograph)):
            print("Warning: Data contains NaN or inf values, will be replaced with 0")
            kymograph_norm = np.nan_to_num(kymograph.copy(), 0)
        else:
            kymograph_norm = kymograph.copy()
        
        # Normalize to uint16 range
        if kymograph_norm.dtype != np.uint16:
            print("Converting to uint16 format...")
            
            # Ensure data is within valid range
            min_val = np.min(kymograph_norm)
            max_val = np.max(kymograph_norm)
            
            if min_val == max_val:
                print("Warning: Data range is 0, will be set to 0 directly")
                kymograph_norm = np.zeros(kymograph_norm.shape, dtype=np.uint16)
            else:
                # Normalize to 0-65535 range
                kymograph_norm = ((kymograph_norm - min_val) * (65535.0 / (max_val - min_val))).astype(np.uint16)
            
            print(f"Converted value range: [{np.min(kymograph_norm)}, {np.max(kymograph_norm)}]")
        
        # Save TIFF file
        print(f"Saving TIFF file: {tiff_path}")
        tifffile.imwrite(tiff_path, kymograph_norm)
        print(f"Successfully saved TIFF file: {tiff_path}")
        
        print(f"\nAll files saved successfully:")
        print(f"- {np_path} (raw data)")
        print(f"- {tiff_path} (16-bit TIFF image)")
        
    except Exception as e:
        print(f"\nError occurred during saving:")
        print(f"Error type: {type(e).__name__}")
        print(f"Error message: {str(e)}")
        print(f"Attempted save path: {base_path}")
        raise

def main():
    try:
        # Set parameters
        start_frame = 0  # Starting frame (0-based)
        num_frames = 100  # Number of frames to process
        line_width = 5  # Vessel line width (pixels)
        num_points = 8 # Number of points user needs to click
        
        # 1. Load TIFF image stack
        filename = r"D:\Data\AES data\processed data\file1plane3_normalized.tif"
        #filename = r"D:\Data\AES data\processed data\mouseB4422_AES_volume_100pt8Hz_130um_1pt8mW_00001_ordered_avg1_MIPs.tif"
        filename = r"E:\AES Data\Deinterleave\file3\mouseL4474_AES_dualPort_50pt4Hz_45-105um_2pt4mW_0pt4mW_00003_plane_6_chunk_13of40.tif"
        #filename = r"D:\Data\AES data\Deinterleave\Mouse1file3\mouse1_AES_dualPort_vol50pt4Hz_3pt6mW_0pt7mW_220-280um_inj2_00003_plane_6_chunk_13of40.tif"
        #filename = r"E:\AES Data\mouseB4422_singlePlane_137Hz_130um_3pt5mW_00001.tif"
        #filename = r"D:\Data\AES data\Deinterleave\Mouse2file1\mouse2_AES_dualPort_vol50pt4Hz_0pt87mW_0pt23mW_55-115um_inj2_00001_plane_1_chunk_13of40.tif"
        filename = r"D:\Data\AES data\processed data\mouseB4422_noAES_dualPort_volume_100pt8Hz_130um_3pt5mW_1pt8mW_00001_ordered_avg1_MIPs.tif"
        filename = r"E:\AES Data\20250604 practice wild type mouse blood vessel volume imaging\mouse imaging\file2\mouse1_AES_dualPort_50pt4Hz_684-732um_3pt4mW_0pt53mW_00002_plane_2_chunk_01of40.tif"
        filename = r"F:\Two photon AES Blood Vessel Data\20251128 bloodvessel volume imaging\mouse_imaging\mouseL2_noAES_370-420um_50pt4Hz_4pt0mW_1pt9mW_00001\MIPs\mouseL2_noAES_370-420um_50pt4Hz_4pt0mW_1pt9mW_00001_MIP_chunk_06of10.tif"
        #filename = r"F:\Two photon AES Blood Vessel Data\20251128 bloodvessel volume imaging\mouse_imaging\mouseL2_noAES_370-420um_sqPix_99pt9Hz_4pt0mW_1pt9mW_00001\MIPs\mouseL2_noAES_370-420um_sqPix_99pt9Hz_4pt0mW_1pt9mW_00001_MIP_chunk_10of20.tif"

        if not os.path.exists(filename):
            print(f"Error: File not found {filename}")
            return
            
        print("Reading image file...")
        image_stack = load_tiff_stack(filename, start_frame, num_frames)
        
        # 2. Display image for user to select vessel path points
        display_img = np.max(image_stack, axis=0)  # Use maximum intensity projection
        
        plt.figure(figsize=(10, 10))
        plt.imshow(display_img, cmap='gray')
        plt.title(f'Click {num_points} points along the vessel', fontname='Arial')
        
        # 3. Get multiple points clicked by user
        print(f"Please click {num_points} points along the vessel")
        points = plt.ginput(num_points, timeout=-1)
        plt.close()
        
        # Convert point format
        vessel_points = np.array(points)
        
        # 4. Interpolate points and display path
        interpolated_points = interpolate_curve(vessel_points)
        
        # Create output directory
        output_dir = "results"
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            
        base_filename = os.path.splitext(os.path.basename(filename))[0]
        save_path = os.path.join(output_dir, 
                               f"{base_filename}_kymograph_curved_frames_{start_frame}-{start_frame+num_frames}_width{line_width}")
        
        # Save and display vessel path
        plot_vessel_path(display_img, vessel_points, interpolated_points, save_path + "_vessel_path.png")
        save_vessel_path(vessel_points, interpolated_points, save_path)
        
        # 5. Generate kymograph
        print("Generating kymograph...")
        kymograph = create_vessel_kymograph(image_stack, vessel_points, line_width=line_width)
        print(f"Kymograph shape: {kymograph.shape}")
        
        # 6. Display results
        plot_kymograph(kymograph, start_frame)
        
        # 7. Save results
        save_kymograph(kymograph, save_path)
        
        # Save visualization plot
        plot_kymograph(kymograph, start_frame, save_path + "_plot.png")
        print(f"Results saved to {output_dir} directory")
        
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        raise

if __name__ == '__main__':
    main() 