import os
import re
import matplotlib.pyplot as plt

experiment = 'experiment4'

def extract_data_from_file(filename):
    with open(filename, 'r') as file:
        content = file.read()
    
    fragments = re.findall(r'-------- Reduction with (\d+) threads --------.*?Elapsed time: ([\d.]+) ms', content, re.DOTALL)
    data = [(int(num_threads), float(elapsed_time)) for num_threads, elapsed_time in fragments]
    return data

def main():
    experiment_results_dir = 'experiment_results'
    experiment_files = [f for f in os.listdir(experiment_results_dir) if f.startswith(experiment)]
    
    results = {}
    for file in experiment_files:
        num_floats = int(re.search(r'_(\d+)_floats\.txt', file).group(1))
        data = extract_data_from_file(os.path.join(experiment_results_dir, file))
        if num_floats not in results:
            results[num_floats] = []
        results[num_floats].extend(data)
    
    for num_floats, data in results.items():
        print(f'Number of floats: {num_floats}')
        for num_threads, elapsed_time in data:
            print(f'  Threads: {num_threads}, Elapsed time: {elapsed_time} ms')

    plt.figure()
    for num_floats, data in results.items():
        threads = [item[0] for item in data]
        times = [item[1] for item in data]
        
        plt.plot(threads, times, marker='o', label=f'{num_floats} floats')
    
    plt.title('Ring: 4 nodes')
    plt.xlabel('Number of Threads')
    plt.ylabel('Elapsed Time (ms)')
    plt.yscale('log')
    plt.legend()
    plt.grid(True)
    plt.xticks([2, 4, 8, 16, 32, 64])
    plt.savefig('img/'+str(experiment)+'.png')
    plt.show()

if __name__ == '__main__':
    main()
    