def print_gc_contents(sequences, gc_func):
    print("\nGC Content:")
    for record in sequences:
        gc = gc_func(record)
        print(record.id, "GC%:", round(gc, 2))
