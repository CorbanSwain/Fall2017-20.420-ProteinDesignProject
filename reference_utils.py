import multiprocessing



# parmap implementation for functions with non-trivial execution time
def funmap(f, q_in, q_out):
    while True:
        i, x = q_in.get()
        if i is None:
            break
        q_out.put((i, f(x)))
    return


def parmap(f, X, nprocs=(multiprocessing.cpu_count() - 1)):
    q_in = multiprocessing.Queue(1)
    q_out = multiprocessing.Queue()
    proc = [multiprocessing.Process(target=funmap, args=(f, q_in, q_out)) for _ in range(nprocs)]
    for p in proc:
        p.daemon = True
        p.start()
    [q_in.put((i, x)) for i, x in enumerate(X)]
    [q_in.put((None, None)) for _ in range(nprocs)]
    res = [q_out.get() for _ in range(len(X))]
    [p.join() for p in proc]
    return [x for i, x in sorted(res)]

def main():
    return 0


if __name__ == "__main__":
    main()

