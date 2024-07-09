import matplotlib.pyplot as plt

def read_lines_from_file(filename):
    lines = []
    with open(filename, 'r') as file:
        for line in file:
            points = list(map(float, line.replace(',', ' ').split()))
            lines.append(points)
    return lines

def parse_points(lines):
    parsed_lines = []
    for points in lines:
        x = points[0::2]
        y = points[1::2]
        parsed_lines.append((x, y))
    return parsed_lines

def plot_lines(lines):
    for line in lines:
        x, y = line
        plt.plot(x, y)

    plt.axis('equal')
    plt.show()

if __name__ == "__main__":
    filename = 'data.txt'
    lines = read_lines_from_file(filename)
    parsed_lines = parse_points(lines)
    plot_lines(parsed_lines)