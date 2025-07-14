# Small file to create a rudimentary progress bar. Generating the dictionary could
# take hours or even days, so its helpful to check in on how it's going.
from sys import stdout
spinner_num = 0
spinner_vals = ["|", "\\", "-", "/"]

def get_spinner_char():
    global spinner_num
    global spinner_vals

    val = spinner_vals[spinner_num]
    spinner_num = (spinner_num + 1) % 4
    return 3 * val
    

def create_pb():
    print() # New Line
    write_pb(0.0)


def refresh_pb(percent):
    clear_pb()
    write_pb(percent)


def clear_pb():
    stdout.write('\r' + " " * 50 + '\r')
    stdout.flush()

def finish_pb():
    clear_pb()
    num_eq = int(20)
    stdout.write("\033[36mProgress: ")
    stdout.write("\033[31m<\033[32m")
    stdout.write(num_eq * "=")
    stdout.write("\033[31m> ")
    stdout.write(f"\033[32mComplete!!!\033[0m\n")
    stdout.flush()


def write_pb(percent):
    num_eq = int(20 * percent)
    stdout.write("\033[36mProgress: ")
    stdout.write("\033[31m<\033[32m")
    stdout.write(num_eq * "=")
    stdout.write((20 - num_eq) * " ")
    stdout.write("\033[31m> ")
    stdout.write(f"\033[36m{100 * percent:2.2f}% ")
    stdout.write("\033[32m" + get_spinner_char() + "\033[0m")
    stdout.flush()


### TESTING ###
if __name__ == "__main__":
    from time import sleep

    create_pb()
    sleep(1)

    for i in range(201):
        refresh_pb(i / 200)
        sleep(0.05)
    
    finish_pb()