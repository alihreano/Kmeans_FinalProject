import os
import subprocess
import numpy as np
import sys


def print_bold(message, end = '\n'):
    sys.stdout.write('\x1b[1;37m' + message.strip() + '\x1b[0m' + end)

def print_fail(message, end='\n'):
    sys.stderr.write('\x1b[1;31m' + message.strip() + '\x1b[0m' + end)


def print_pass(message, end='\n'):
    sys.stdout.write('\x1b[1;32m' + message.strip() + '\x1b[0m' + end)


def compare_terminal_output_to_file(test_name, command, file_path):
    p1 = subprocess.Popen(command.split(), stdout=subprocess.PIPE)
    p2 = subprocess.Popen(["diff", file_path, "-"],
                          stdin=p1.stdout, stdout=subprocess.PIPE)
    output, _ = p2.communicate()
    if len(output) == 0:
        print_pass(f"{test_name}\tpass")
    else:
        print_fail(f"{test_name}\tfail")


def test_wam_ddg_gl_C():
    for func in ["wam", "ddg", "gl"]:
        for i in range(1, 11):
            test_name = f"{func} test {i}"
            command = f"./spkmeans {func} test/test_batch/test{i}.txt"
            expected = f"test/ans_batch/test{i}_{func}.txt"
            compare_terminal_output_to_file(test_name, command, expected)


def test_jacobi_C():
    for i in range(1, 11):
        test_name = f"jacobi test {i}"
        command = f"./spkmeans jacobi test/test_batch/test{i}_j.txt"
        expected = f"test/ans_batch/test{i}_j_ans.txt"
        compare_terminal_output_to_file(test_name, command, expected)


def test_wam_ddg_gl_python():
    for func in ["wam", "ddg", "gl"]:
        for i in range(1, 11):
            test_name = f"{func} test {i}"
            command = f"python3 spkmeans.py {func} test/test_batch/test{i}.txt"
            expected = f"test/ans_batch/test{i}_{func}.txt"
            compare_terminal_output_to_file(test_name, command, expected)


def test_jacobi_python():
    for i in range(1, 11):
        test_name = f"jacobi test {i}"
        command = f"python3 spkmeans.py jacobi test/test_batch/test{i}_j.txt"
        expected = f"test/ans_batch/test{i}_j_ans.txt"
        compare_terminal_output_to_file(test_name, command, expected)


def test_spk_python():
    for i in range(1, 11):
        test_name = f"spk test {i}"
        command = f"python3 spkmeans.py spk test/test_batch/test{i}.txt"
        expected = f"test/ans_batch/test{i}_spk.txt"
        compare_terminal_output_to_file(test_name, command, expected)


def filter_valgrind_output(output):
    def filter_function(line): return line if "in use at exit" in line else ""
    filtered = filter(filter_function, output.splitlines())
    return list(filtered)[0].split("     ")[-1]


def test_leaks():
    for func in ["wam", "ddg", "gl"]:
        for i in range(1, 11):
            valgrind_command = ["valgrind", "--leak-check=yes",
                                "./spkmeans", func, f"test/test_batch/test{i}.txt"]
            output = subprocess.check_output(
                valgrind_command, text=True, stderr=subprocess.STDOUT)
            relevant_output = filter_valgrind_output(output)
            passed = relevant_output == 'in use at exit: 0 bytes in 0 blocks'
            out = f"{func} test {i}\t" + relevant_output
            print_pass(out) if passed else print_fail(out)

    for i in range(1, 11):
        valgrind_command = ["valgrind", "--leak-check=yes",
                            "./spkmeans", "jacobi", f"test/test_batch/test{i}_j.txt"]
        output = subprocess.check_output(
            valgrind_command, text=True, stderr=subprocess.STDOUT)
        relevant_output = filter_valgrind_output(output)
        passed = relevant_output == 'in use at exit: 0 bytes in 0 blocks'
        out = f"jacobi test {i}\t" + relevant_output
        print_pass(out) if passed else print_fail(out)


def test_all():
    os.system('make')
    os.system('clear')
    print_bold(f"{'-'*3}Testing Python{'-'*3}")
    test_wam_ddg_gl_python()
    test_jacobi_python()
    test_spk_python()
    print_bold(f"{'-'*3}Testing C{'-'*3}")
    test_wam_ddg_gl_C()
    test_jacobi_C()
    print_bold(f"{'-'*15}Testing Memory Leaks:{'-'*15}")
    test_leaks()


if __name__ == "__main__":
    test_all()

        
