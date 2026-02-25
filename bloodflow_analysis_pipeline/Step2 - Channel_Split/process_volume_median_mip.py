import numpy as np
import tifffile
from pathlib import Path
import argparse
from tqdm import tqdm
import multiprocessing as mp
from math import ceil
import os
from scipy.ndimage import median_filter

def load_plane_chunk(plane_file_path):
    """
    Load a single plane chunk file.
    
    Args:
        plane_file_path (str): Path to the plane chunk file
        
    Returns:
        numpy.ndarray: Loaded image data
    """
    try:
        with tifffile.TiffFile(plane_file_path) as tif:
            # 检查是否有多个页面
            if len(tif.pages) > 1:
                # 手动读取所有页面
                all_pages = []
                for i in range(len(tif.pages)):
                    page_data = tif.pages[i].asarray()
                    all_pages.append(page_data)
                data = np.stack(all_pages, axis=0)
            else:
                # 单页文件，直接读取
                data = tif.asarray()
            
            print(f"Loaded {plane_file_path}: shape={data.shape}, pages={len(tif.pages)}")
            return data
    except Exception as e:
        print(f"Error loading {plane_file_path}: {e}")
        return None

def process_chunk_mip(args):
    """
    Process a single chunk to create MIP across all 6 planes.
    
    Args:
        args: tuple containing (input_dir, base_name, chunk_id, total_chunks, output_dir)
    """
    input_dir, base_name, chunk_id, total_chunks, output_dir = args
    
    # Load all 6 planes for this chunk
    plane_data = []
    for plane in range(1, 7):  # 6 planes
        plane_file = input_dir / f"{base_name}_plane_{plane}_chunk_{chunk_id+1:02d}of{total_chunks:02d}.tif"
        
        if not plane_file.exists():
            print(f"Warning: {plane_file} not found, skipping chunk {chunk_id+1}")
            return None
        
        data = load_plane_chunk(str(plane_file))
        if data is None:
            return None
        
        # Apply median filter to each frame in the plane
        # The filter is 2D, so we specify a size of (1, 3, 3) to operate on (height, width)
        # for each frame, leaving the frame axis (time) untouched.
        filtered_data = median_filter(data, size=(1, 3, 3))
        plane_data.append(filtered_data)
    
    # Stack all planes and calculate MIP
    # plane_data is a list of 6 arrays, each with shape (frames, height, width)
    stacked_volumes = np.stack(plane_data, axis=1)  # Shape: (frames, 6_planes, height, width)
    
    # Calculate MIP across all planes for each frame
    mips = np.max(stacked_volumes, axis=1)  # Shape: (frames, height, width)
    
    # Save the MIP chunk
    output_file = output_dir / f"{base_name}_MIP_chunk_{chunk_id+1:02d}of{total_chunks:02d}.tif"
    tifffile.imwrite(str(output_file), mips, compression='zlib')
    
    return output_file

def create_mip_video(input_dir, output_dir, base_name, num_chunks=40):
    """
    Create maximum intensity projection (MIP) video from deinterleaved plane data.
    
    Args:
        input_dir (str): Directory containing the deinterleaved plane files
        output_dir (str): Directory to save the MIP video files
        base_name (str): Base name of the files (without extension)
        num_chunks (int): Number of chunks to process (default: 40)
    """
    input_dir = Path(input_dir)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    
    print(f"Processing MIP from: {input_dir}")
    print(f"Output directory: {output_dir}")
    print(f"Base name: {base_name}")
    print(f"Number of chunks: {num_chunks}")
    
    # Prepare arguments for parallel processing
    chunk_args = []
    for chunk_id in range(num_chunks):
        chunk_args.append((input_dir, base_name, chunk_id, num_chunks, output_dir))
    
    # Process chunks in parallel
    print("\nProcessing chunks...")
    with mp.Pool(min(8, num_chunks)) as pool:
        results = list(tqdm(
            pool.imap(process_chunk_mip, chunk_args),
            total=num_chunks,
            desc="Processing MIP chunks"
        ))
    
    # Filter out None results (failed chunks)
    successful_results = [r for r in results if r is not None]
    
    print(f"\nSuccessfully processed {len(successful_results)} chunks")
    print(f"MIP files saved to: {output_dir}")
    
    return successful_results

def main():
    parser = argparse.ArgumentParser(
        description='Create maximum intensity projection (MIP) video from deinterleaved plane data'
    )
    parser.add_argument('input_dir', help='Directory containing the deinterleaved plane files')
    parser.add_argument('output_dir', help='Directory to save the MIP video files')
    parser.add_argument('--base-name', help='Base name of the files (if not specified, will be inferred from input directory)')
    parser.add_argument('--num-chunks', type=int, default=40,
                      help='Number of chunks to process (default: 40)')
    
    args = parser.parse_args()
    
    # Determine base name if not provided
    if args.base_name is None:
        # Try to infer from the first plane file
        input_dir = Path(args.input_dir)
        plane_files = list(input_dir.glob("*_plane_1_chunk_01of*.tif"))
        if plane_files:
            filename = plane_files[0].stem
            # Extract base name by removing "_plane_1_chunk_01ofXX"
            base_name = filename.replace("_plane_1_chunk_01of", "").rsplit("_", 1)[0]
        else:
            print("Error: Could not find plane files to infer base name. Please specify --base-name.")
            return
    else:
        base_name = args.base_name
    
    print(f"Using base name: {base_name}")
    
    # Create MIP chunks
    mip_files = create_mip_video(args.input_dir, args.output_dir, base_name, args.num_chunks)
    
    print("\nProcessing complete!")

if __name__ == "__main__":
    # Example usage
    input_dir = r"F:\Two photon AES Blood Vessel Data\20251128 bloodvessel volume imaging\mouse_imaging\mouseL2_noAES_370-420um_sqPix_99pt9Hz_4pt0mW_1pt9mW_00001"
    output_dir = r"F:\Two photon AES Blood Vessel Data\20251128 bloodvessel volume imaging\mouse_imaging\mouseL2_noAES_370-420um_sqPix_99pt9Hz_4pt0mW_1pt9mW_00001\MIPs"
    base_name = "mouseL2_noAES_370-420um_sqPix_99pt9Hz_4pt0mW_1pt9mW_00001"
    
    if Path(input_dir).exists():
        print("Running example...")
        # Create MIP chunks
        mip_files = create_mip_video(input_dir, output_dir, base_name, num_chunks=20)
    else:
        print(f"Example input directory not found: {input_dir}")
        print("Please run with command line arguments:")
        print("python process_volume_mip.py <input_dir> <output_dir> [--base-name <name>]") 