import tifffile
import numpy as np
from pathlib import Path
import argparse
from tqdm import tqdm
import multiprocessing as mp
from math import ceil
from tqdm.auto import tqdm as auto_tqdm
import os

def shift_image_up(image, pixels=2):
    """
    Shift image up by specified number of pixels.
    The bottom rows will be filled with the last row's values.
    
    Args:
        image: 2D numpy array
        pixels: number of pixels to shift up
    """
    if pixels == 0:
        return image
    
    # Create shifted image
    shifted = np.roll(image, -pixels, axis=0)
    # Fill the bottom rows with the last valid row's values
    shifted[-pixels:] = shifted[-pixels-1]
    return shifted

def print_frame_indices(num_cycles=3):
    """
    Print the frame indices for each plane in the first few cycles.
    """
    # Arrays indicate which frame goes to which plane
    # e.g., frame 2 goes to plane 1, frame 4 goes to plane 2, etc.
    forward_frames = [1,3,5,2,4,6]    # 1-based indexing
    backward_frames = [11,9,7,12,10,8] # 1-based indexing
    
    print("\nFrame mapping for first few cycles (1-based indexing):")
    for plane in range(1, 7):
        frames = []
        for cycle in range(num_cycles):
            cycle_start = cycle * 12  # 1-based indexing
            # Forward scan frame for this plane
            forward_frame = cycle_start + forward_frames[plane-1]
            frames.append(forward_frame)
            # Backward scan frame for this plane
            backward_frame = cycle_start + backward_frames[plane-1]
            frames.append(backward_frame)
        print(f"Plane {plane} comes from frames: {frames}")

def process_cycles(args):
    """
    Process a chunk of scan cycles.
    
    Args:
        args: tuple containing (input_path, output_dir, start_cycle, end_cycle, y_offset, chunk_id, total_chunks)
    """
    input_path, output_dir, start_cycle, end_cycle, y_offset, chunk_id, total_chunks = args
    
    # Arrays indicate which frame goes to which plane
    forward_frames = [2,4,6,1,3,5]    # 1-based indexing
    backward_frames = [12,10,8,11,9,7] # 1-based indexing
    
    # Convert to 0-based indexing for Python
    forward_frames = [x-1 for x in forward_frames]
    backward_frames = [x-1 for x in backward_frames]
    
    # Create output directory if it doesn't exist
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Get the base filename without extension
    base_name = Path(input_path).stem
    
    # Create output files for this chunk (one for each plane)
    output_files = []
    for i in range(6):  # 6 planes
        output_path = output_dir / f"{base_name}_plane_{i+1}_chunk_{chunk_id+1:02d}of{total_chunks:02d}.tif"
        output_files.append(tifffile.TiffWriter(str(output_path)))
    
    try:
        # Open the TIFF file
        with tifffile.TiffFile(input_path) as tif:
            # Create progress bar for this chunk
            pbar_desc = f"Chunk {chunk_id+1:2d}/{total_chunks}"
            with auto_tqdm(total=end_cycle - start_cycle,
                          desc=pbar_desc,
                          position=chunk_id % 8,  # Show 8 progress bars at a time
                          leave=True) as pbar:
                
                # Process assigned cycles
                for cycle in range(start_cycle, end_cycle):
                    cycle_start = cycle * 12  # 12 frames per cycle
                    
                    # Process forward scan frames
                    for plane_idx in range(6):
                        frame_idx = forward_frames[plane_idx]
                        page_idx = cycle_start + frame_idx
                        
                        # Read and write the page
                        page_data = tif.pages[page_idx].asarray()
                        output_files[plane_idx].write(
                            page_data,
                            photometric='minisblack',
                            compression='zlib'
                        )
                    
                    # Process backward scan frames
                    for plane_idx in range(6):
                        frame_idx = backward_frames[plane_idx]
                        page_idx = cycle_start + frame_idx
                        
                        # Read the page, flip Y direction, and apply Y offset for backward scan
                        page_data = tif.pages[page_idx].asarray()
                        page_data = np.flipud(page_data)  # Flip Y direction
                        page_data = shift_image_up(page_data, y_offset)  # Apply Y offset
                        
                        # Write to the output file
                        output_files[plane_idx].write(
                            page_data,
                            photometric='minisblack',
                            compression='zlib'
                        )
                    
                    pbar.update(1)
    finally:
        # Close all output files
        for f in output_files:
            f.close()

def deinterleave_bidirectional(input_path, output_dir, y_offset=2, num_chunks=40):
    """
    Deinterleave a large multi-page TIFF file from bidirectional scanning using parallel processing.
    Data will be split into num_chunks parts and saved separately.
    
    Args:
        input_path (str): Path to the input TIFF file
        output_dir (str): Directory to save the deinterleaved files
        y_offset (int): Number of pixels to shift up for backward scan compensation
        num_chunks (int): Number of chunks to split the data into
    """
    print_frame_indices()
    
    # Get total number of cycles
    with tifffile.TiffFile(input_path) as tif:
        total_pages = len(tif.pages)
        frames_per_cycle = 12
        total_cycles = total_pages // frames_per_cycle
    
    print(f"\nProcessing file: {input_path}")
    print(f"Total cycles: {total_cycles}")
    print(f"Splitting into {num_chunks} chunks")
    
    # Calculate cycles per chunk
    base_cycles_per_chunk = total_cycles // num_chunks
    extra_cycles = total_cycles % num_chunks
    
    # Prepare chunk arguments
    chunk_args = []
    start_cycle = 0
    
    for i in range(num_chunks):
        # Add an extra cycle to some chunks if total_cycles is not evenly divisible
        cycles_for_this_chunk = base_cycles_per_chunk + (1 if i < extra_cycles else 0)
        end_cycle = start_cycle + cycles_for_this_chunk
        
        if cycles_for_this_chunk > 0:  # Only create work for non-empty chunks
            chunk_args.append((input_path, output_dir, start_cycle, end_cycle, y_offset, i, num_chunks))
        
        start_cycle = end_cycle
    
    # Clear some space for progress bars (showing 8 at a time)
    print("\n" * 10)
    
    # Create a pool of workers and process the chunks in parallel
    with mp.Pool(min(8, num_chunks)) as pool:  # Use at most 8 processes at a time
        pool.map(process_cycles, chunk_args)
    
    print("\nDeinterleaving complete!")
    print(f"Output files have been saved to: {output_dir}")
    base_name = Path(input_path).stem
    print(f"Created {6 * num_chunks} files (6 planes × {num_chunks} chunks)")
    print("\nFile naming pattern:")
    print(f"  {base_name}_plane_X_chunk_YYofZZ.tif")
    print("where:")
    print("  X: plane number (1-6)")
    print(f"  YY: chunk number (01-{num_chunks:02d})")
    print(f"  ZZ: total chunks ({num_chunks:02d})")

def main():
    parser = argparse.ArgumentParser(
        description='Deinterleave a multi-page TIFF file from bidirectional scanning'
    )
    parser.add_argument('input_path', help='Path to the input TIFF file')
    parser.add_argument('output_dir', help='Directory to save the deinterleaved files')
    parser.add_argument('--y-offset', type=int, default=3,
                      help='Number of pixels to shift up for backward scan compensation (default: 2)')
    parser.add_argument('--num-chunks', type=int, default=40,
                      help='Number of chunks to split the data into (default: 40)')
    
    args = parser.parse_args()
    deinterleave_bidirectional(args.input_path, args.output_dir, 
                             args.y_offset, args.num_chunks)

if __name__ == '__main__':
    main() 