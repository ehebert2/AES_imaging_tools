import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import splprep, splev
from skimage.measure import profile_line
import tifffile
import os
import re
from pathlib import Path

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

def load_segments_from_txt(filepath):
    """
    Load coordinates of multiple vessel segments from a txt file.
    File format should be one point per line, with x and y coordinates separated by space or comma.
    Segments are separated by lines starting with '#'.
    Example:
    # Segment 1
    100,200
    150,250
    # Segment 2
    200,300
    210,310
    
    Parameters:
    filepath : str
        Path to the txt file containing coordinates
        
    Returns:
    list of numpy.ndarray
        A list where each element is a numpy array of shape (N, 2), representing coordinate points of a vessel segment.
    """
    segments = []
    current_segment = []
    try:
        with open(filepath, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    if current_segment:
                        segments.append(np.array(current_segment))
                        current_segment = []
                elif line:
                    try:
                        coords = [float(x) for x in line.replace(',', ' ').split()]
                        if len(coords) == 2:
                            # Swap X and Y coordinates
                            current_segment.append([coords[1], coords[0]])
                        else:
                            print(f"Warning: Ignoring incorrectly formatted line: '{line}'")
                    except ValueError:
                        print(f"Warning: Ignoring unparseable line: '{line}'")
            
            if current_segment:
                segments.append(np.array(current_segment))
        
        print(f"Loaded {len(segments)} vessel segment(s) from {filepath}")
        if not segments:
            print("Warning: Failed to parse any valid vessel segments from file. Please check file format.")
        else:
            for i, seg in enumerate(segments):
                print(f"  - Segment {i+1}: {len(seg)} points")
            
        return segments
    except Exception as e:
        print(f"Error loading point file: {e}")
        raise

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

    # Build a (n,2) array of points from x,y
    pts = np.column_stack([x, y])
    
    # 1) Remove exact duplicate points
    # (keep order and first occurrence)
    _, idx = np.unique(pts, axis=0, return_index=True)
    pts_unique = pts[np.sort(idx)]
    
    # 2) Ensure no NaN/Inf
    mask = np.isfinite(pts_unique).all(axis=1)
    pts_unique = pts_unique[mask]
    
    n = len(pts_unique)
    if n < 2:
        raise ValueError("Need at least 2 unique, finite points for a curve.")
    
    # 3) Choose a safe k given unique points
    # cubic (k=3) needs at least 4 unique points; otherwise fallback
    k = min(3, max(1, n - 1))
    # Use splprep for spline interpolation
    tck, u = splprep([pts_unique[:,0], pts_unique[:,1]], s=0, k=k)


    # # Use splprep for spline interpolation
    # tck, u = splprep([x, y], s=0, k=min(3, len(points)-1))
    
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
                                 mode='constant',
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

# def plot_vessel_path_old(image, all_points, all_interpolated_points=None, save_path=None):
#     """
#     Display vessel path on image with optional saving
    
#     Parameters:
#     image : numpy.ndarray
#         Image to display
#     all_points : list of numpy.ndarray
#         List containing original points of all segments
#     all_interpolated_points : list of numpy.ndarray, optional
#         List containing interpolated points of all segments
#     save_path : str, optional
#         Save path (without extension)
#     """
#     plt.figure(figsize=(10, 10))
#     plt.imshow(image, cmap='gray')
    
#     # Plot all segments
#     for i, points in enumerate(all_points):
#         # Plot user-selected points
#         plt.plot(points[:, 0], points[:, 1], 'o', label=f'Segment {i+1} points')
    
#     if all_interpolated_points:
#         for i, interpolated_points in enumerate(all_interpolated_points):
#             # If interpolated points exist, plot them as well
#             plt.plot(interpolated_points[:, 0], interpolated_points[:, 1], '-', 
#                      label=f'Segment {i+1} path', linewidth=1)
    
#     plt.legend()
#     plt.title('Vessel Paths from File')
    
#     if save_path:
#         plt.savefig(save_path, dpi=300, bbox_inches='tight')
#     plt.show()

def plot_vessel_path(image, all_points, all_interpolated_points=None, save_path=None):
    plt.figure(figsize=(10, 10))
    plt.imshow(image, cmap='gray')

    # Generate a color array with as many colors as segments
    num_segments = len(all_points)
    colors = plt.cm.hsv(np.linspace(0, 1, num_segments))  # hsv gives a rainbow-like spread
    colors2 = plt.cm.viridis(np.linspace(0, 1, num_segments))

    for i, points in enumerate(all_points):
        color = colors[i]
        plt.plot(points[:, 0], points[:, 1], 'o', color=color, label=f'Segment {i+1} points')

    if all_interpolated_points:
        for i, interpolated_points in enumerate(all_interpolated_points):
            color = colors2[(41*i)%num_segments]
            plt.plot(interpolated_points[:, 0], interpolated_points[:, 1], '-', 
                     color=color, label=f'Segment {i+1} path', linewidth=1)

    plt.legend()
    plt.title('Vessel Paths from File')

    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
    plt.show()


# def save_vessel_path(points, interpolated_points, base_path):
#     """
#     Save vessel path data
    
#     Parameters:
#     points : numpy.ndarray
#         User-selected points
#     interpolated_points : numpy.ndarray
#         Interpolated points
#     base_path : str
#         Base save path (without extension)
#     """
#     # Save original points and interpolated points
#     path_data = {
#         'original_points': points,
#         'interpolated_points': interpolated_points
#     }
#     np.savez(f"{base_path}_path.npz", **path_data)
#     print(f"- {base_path}_path.npz (vessel path data)")
    


def save_vessel_path(points, interpolated_points, base_path):
    """
    Save vessel path data.

    Parameters:
    points : numpy.ndarray
        User-selected points
    interpolated_points : numpy.ndarray
        Interpolated points
    base_path : str or pathlib.Path
        Base save path (without extension)
    """
    # Normalize to Path (handles strings with backslashes or slashes)
    base_path = Path(base_path)

    # Build the target .npz path: base + "_path.npz"
    npz_path = base_path.with_name(base_path.name + "_path").with_suffix(".npz")

    # Ensure parent directory exists
    npz_path.parent.mkdir(parents=True, exist_ok=True)

    # Save data
    path_data = {
        "original_points": points,
        "interpolated_points": interpolated_points,
    }
    np.savez(npz_path, **path_data)

    print(f"- Saved vessel path: {npz_path.resolve()}")



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
        # --- Parameter Settings ---
        
        # 1. TIF file path
        # Please modify this path to your TIF file path
          
        # tif_filename = r"E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\processed_videos\MAX_mouseA2_subY_aM0pt8_50pt4Hz_890-940um_flip_video4_vol9601-19200.tif"
        # tif_filename = r"E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\processed_videos\MAX_mouseL4_video4_processed_red_vol1-23000_motCorr.tif"
        tif_filename = r"E:\Volumetric blood vessel imaging data\20251220 bloodvessel volume imaging\processed_videos\MAX_mouseL4_video4_processed_red_vol23001-46000_motCorr.tif"
        
        # 2. Coordinate file path
        # Please modify this path to your coordinate file path
        points_filename = "SegmentsPrompt.txt" # Example filename, please ensure this file exists

        # x_FoV_shift = 1   # 200 Hz video sq pix
        # y_FoV_shift = 169

        # x_FoV_shift = 1   # original video
        # y_FoV_shift = 116

        # x_FoV_shift = 1    # 200 Hz video 1.8x subYsample
        # y_FoV_shift = 136
        # x_FoV_shift = 2    # 100 Hz video 1.96x subYsample
        # y_FoV_shift = 45
        # x_FoV_shift = 1    # 150 Hz video 1.7x subYsample
        # y_FoV_shift = 114
        # x_FoV_shift = 1    # 150 Hz video 1.7x subYsample, 20251205 700 um
        # y_FoV_shift = 120
        # x_FoV_shift = -13    # 150 Hz video 1.7x subYsample, 20251220 400 um
        # y_FoV_shift = 120
        x_FoV_shift = -5    # 150 Hz video 1.7x subYsample, 20251220 400 um video4
        y_FoV_shift = 117
        # x_FoV_shift = -5+6   # 150 Hz video 1.7x subYsample, 20251220 400 um video4 with motion
        # y_FoV_shift = 117+4
        # x_FoV_shift = -5+4   # 150 Hz video 1.7x subYsample, 20251220 400 um video4 with motion
        # y_FoV_shift = 117+0

        # 3. Kymograph parameters
        start_frame = 16000  # Starting frame (0-based)
        num_frames = 7000  # Number of frames to process
        line_width = 5    # Vessel line width (pixels, must be odd)
        
        # --- Main Program ---

        # Check if files exist
        if not os.path.exists(tif_filename):
            print(f"Error: TIF file not found {tif_filename}")
            return
        if not os.path.exists(points_filename):
            print(f"Error: Coordinate file not found {points_filename}")
            print("Please create a txt file containing vessel path coordinates, with format 'x y' per line, segments separated by blank lines")
            return
            
        # 1. Load TIFF image stack
        print("Reading image file...")
        image_stack = load_tiff_stack(tif_filename, start_frame, num_frames)
        
        # 2. Load multiple vessel segment path points from file
        vessel_segments = load_segments_from_txt(points_filename)
        
        if not vessel_segments:
            print("Error: No vessel segments loaded from coordinate file.")
            return

        # Create output directory (Cannot exceed 260 characters!)
        # output_dir = "results"
        output_dir = r"E:\Volumetric blood vessel imaging data\testFolder11"
        os.makedirs(output_dir, exist_ok=True)  # Changed from if/else to always create with exist_ok
            
        base_filename = os.path.splitext(os.path.basename(tif_filename))[0]

        all_original_points_shifted = []
        all_interpolated_points = []

        # 3. Loop through each vessel segment
        for i, segment_points in enumerate(vessel_segments):
            print(f"\n--- Processing segment {i+1}/{len(vessel_segments)} ---")

            # manually correct for FoV_shift
            y_FoV_shift_adj = y_FoV_shift
            x_FoV_shift_adj = x_FoV_shift
            # if i== 11+12-1:
            #     x_FoV_shift_adj = x_FoV_shift - 2
            # elif i == 12+12-1:
            #     x_FoV_shift_adj = x_FoV_shift - 3
            # elif i == 24+1-1:
            #     x_FoV_shift_adj = x_FoV_shift + 3
                
            # Apply offset to coordinates
            segment_points_shifted = segment_points - np.array([x_FoV_shift_adj, y_FoV_shift_adj])
            
            # Interpolate points
            interpolated_points = interpolate_curve(segment_points_shifted)
            
            all_original_points_shifted.append(segment_points_shifted)
            all_interpolated_points.append(interpolated_points)

            # Create unique save path prefix for current segment
            save_path_prefix = os.path.join(output_dir, 
                                   f"{base_filename}_kymo_segment_{i+1}_frames_{start_frame}-{start_frame+num_frames}_width{line_width}")
            
            # print("DEBUG save prefix:", save_path_prefix)
            
            # prefix = Path(save_path_prefix)
            # print("About to save to:", prefix)
            # print("Parent exists?:", prefix.parent.exists())
            
            # # Optional: create parent folder here too (belt-and-suspenders)
            # prefix.parent.mkdir(parents=True, exist_ok=True)
            

            # npz_path = Path(f"{save_path_prefix}_path.npz").resolve()
            # print("Absolute target path:", npz_path)
            # print("Path length:", len(str(npz_path)))


            # Save path data
            save_vessel_path(segment_points_shifted, interpolated_points, save_path_prefix)
            
            # Generate kymograph
            print("Generating kymograph...")
            kymograph = create_vessel_kymograph(image_stack, segment_points_shifted, line_width=line_width)
            print(f"Kymograph shape: {kymograph.shape}")
            
            # Display and save results
            plot_kymograph(kymograph, start_frame, save_path_prefix + "_plot.png")
            save_kymograph(kymograph, save_path_prefix)

        # Display all loaded paths on one image for verification
        display_img = np.max(image_stack, axis=0) # Use maximum intensity projection
        combined_path_save_path = os.path.join(output_dir, f"{base_filename}_kymograph_from_file_ALL_SEGMENTS_vessel_path.png")
        plot_vessel_path(display_img, all_original_points_shifted, all_interpolated_points, combined_path_save_path)
        
        print(f"\nProcessing complete! Results for all {len(vessel_segments)} vessel segment(s) have been saved to {output_dir} directory")
        
    except Exception as e:
        print(f"\nError occurred during program execution: {str(e)}")
        # raise # Uncomment this line if you need to see detailed error stack trace

if __name__ == '__main__':
    main() 