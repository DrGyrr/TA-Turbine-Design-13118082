
import inspect
import sys
print("this print is outside of function")
def purinto():
    if inspect.stack()[1].function == 'puri':
        print("mod1test dipakai di mod2test.puri")
    else:
        print("nope")
        print(inspect.stack()[1].function)