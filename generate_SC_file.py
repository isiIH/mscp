import random

def generate_set_cover_file(filename, n=500, m=200, num_groups=2):
    universe = list(range(1, n + 1))
    subsets = [[] for _ in range(m)]

    #divide universe in groups
    group_size = n // num_groups
    groups = [universe[i * group_size: (i + 1) * group_size] for i in range(num_groups)]

    left_subsets = m
    for i, group in enumerate(groups):
        n_subsets = random.randint(1, min(left_subsets,  m // num_groups))
        if(i == len(groups) - 1):
            n_subsets = left_subsets

        # Add every element inside a random subset
        for e in group:
            rand_subset = random.randint(m - left_subsets, (m - left_subsets) + n_subsets - 1)
            subsets[rand_subset].append(e)

        # Add more elements to random subsets
        for j,subset in enumerate(subsets[m - left_subsets:(m - left_subsets) + n_subsets]):
            subset_size = random.randint(1, len(group) // 4)
            subset_size -= len(subset)
            if subset_size > 0:
                subset += random.sample(group, subset_size)
            # print(f"Group {i} ({j}): {len(subset)} - {subset}")

        left_subsets -= n_subsets

    with open(filename, "w") as file:
        file.write(f"{n} {m}\n")
        for subset in subsets:
            file.write("1 1 " + " ".join(map(str, subset)) + "\n")

    print(f"Dataset guardado en {filename}")

generate_set_cover_file("test/ex9.txt", n=10000, m=20000, num_groups=16)