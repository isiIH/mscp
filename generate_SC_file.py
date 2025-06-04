import random

def generate_set_cover_file(filename, n=500, m=200, num_groups=2):
    universe = list(range(1, n + 1))
    subsets = [[] for _ in range(m)]
    groups = [[] for _ in range(num_groups)]

    #divide universe in groups
    for i in range(num_groups):
        groups[i].append(universe[i])
    # Randomly assign the rest of the elements to groups
    for e in universe[num_groups:]:
        group = random.randint(0, num_groups - 1)
        groups[group].append(e)

    left_subsets = m
    for i, group in enumerate(groups):
        n_subsets = random.randint(1, min(left_subsets,  m // num_groups))
        if(i == len(groups) - 1):
            n_subsets = left_subsets
        # print(f"Group {i}: {len(group)} elements, {n_subsets} subsets")

        # Add every element inside a random subset
        for e in group:
            element_subsets = random.randint(1, n_subsets)
            # Randomly select a list of subsets to add the element to
            ss = random.sample(subsets[m - left_subsets:(m - left_subsets) + n_subsets], element_subsets)

            [s.append(e) for s in ss]
   
        for j,subset in enumerate(subsets[m - left_subsets:(m - left_subsets) + n_subsets]):
            if len(subset) == 0:
                subset += random.sample(group, random.randint(1, len(group) // 2 + 1))
            subset.sort()
            # print(f"({j}): {len(subset)} - {subset}")

        left_subsets -= n_subsets

    assert(sum([len(group) for group in groups]) == n)

    with open(filename, "w") as file:
        file.write(f"{n} {m}\n")
        for subset in subsets:
            file.write("1 1 " + " ".join(map(str, subset)) + "\n")

    print(f"Dataset guardado en {filename}")

generate_set_cover_file("test/test01.txt", n=5000, m=10000, num_groups=4)