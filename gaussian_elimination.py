#!/bin/python3


# Malte Lau Petersen
# maltelau@protonmail.com
# 18-01-2018



### Helper function

def first_nonzero(li):
    # https://stackoverflow.com/questions/19502378/python-find-first-instance-of-non-zero-number-in-list
    return(next((i for i, x in enumerate(li) if x), None))

### Class definition

class equation_system():
    # solver for linear systems of equations
    # how to use:
    # 1) instantiate:
    # >> x = equation_system([[1,3,2],  # left hand side
    #                        [4,0,1],
    #                        [2,1,1]],
    #                        [14,       # right hand side
    #                          5,
    #                          0]) 
    # 2) chose an algorithm:
    # >> x.to_trappeform()
    # >> x.to_reduceret_trappeform()

    def __init__(self, lhs, rhs, debug = False, print_log = True, full_log = False, circle_pivots = False):
        # lhs: matrix
        # rhs: list

        self.lhs = lhs
        self.rhs = rhs
        self.debug = debug
        self.print_log = print_log
        self.full_log = full_log
        self.leading = None

        if circle_pivots:
            self.mark_leading = "({})"
        else:
            self.mark_leading = "{}"

        self.columns = len(lhs)
        self.rows = len(rhs)
        
        self.log = self.__str__() + "\n"
        
    def __str__(self):
        self.calculate_leading()
        # print the matrix
        out = "\n"
        for row in range(self.rows):
            out += "|"
            for col in range(self.columns):
                if col == self.leading[row]:
                    out += "{:>6.5}, ".format(self.mark_leading.format(self.lhs[row][col]))
                elif type(self.lhs[row][col]) is int:
                    out += "{:>6}, ".format(self.lhs[row][col])
                else:
                    out += "{:>6.3}, ".format(float(self.lhs[row][col]))
            out = out[:-2]
            if type(self.rhs[row]) is int:
                out +=" |{:>6} |\n".format(self.rhs[row])
            else:
                out +=" |{:>6.3} |\n".format(self.rhs[row])
        out = out[:-1]
        return(out)

    def to_triangular(self):
        if self.debug: print("\n============== Beginning algorithm to get triangular  ===========") 
        working_on_row = 0
        while True:
            self.calculate_leading()
            if self.debug: print(self)
            
            # step 1
            first_zeroes = [first_nonzero(self.lhs[row]) for row in range(working_on_row, self.rows)]
            pivot_column = min(x for x in first_zeroes if x is not None)
            # step 2
            if self.lhs[working_on_row][pivot_column] == 0:
                self.swap_rows(working_on_row, first_zeroes.index(pivot_column) + working_on_row)

                
            # step 3: brug B E G
            for row in range(working_on_row+1, self.rows):
                # find all rows below the pivot that are not zero
                if self.lhs[row][pivot_column] == 0: continue
                # else: what's B in row2 = row2 + B * row1?
                self.multiply_rows(working_on_row, row, pivot_column)
                
            # step 4:
            working_on_row += 1

            if self.is_triangular():
                self.log += "\nNow in triangular form.\n"
                self.log += self.__str__() + "\n"
                if self.print_log:
                    print("\n=========================== LOG =================================")
                    print(self.log)
                return self
        
    def is_triangular(self):
        # check if the system is currently in triangular
        row_is_null = []
        for row in range(self.rows):
            row_is_null.append(all(self.lhs[row][col] == 0 for col in range(self.columns)))
        
        self.calculate_leading()

            
        # condition 1
        if row_is_null != sorted(row_is_null):
            # There's a nonzero row above a zero row
            if self.debug: print("There's a nonzero row below a zero row.")
            self.swap_rows(row_is_null.index(True), self.rows-1)
            return False

        for row in range(self.rows):
            if self.leading[row] is None:
                # no leading coefficient in this row
                continue
            
            # condition 2
            if not all(self.lhs[x][self.leading[row]] == 0 for x in range(row+1, self.rows)):
                if self.debug: print("Not all coefficients below the leading coefficient [{}, {}] are zero.".format(row+1, self.leading[row]+1))
                return False
            
            # condition 3
            if not all(self.leading[row] > self.leading[x] for x in range(row)):
                if self.debug:
                    print("The leading coefficient at [{}, {}] is not to the right of all leading coefficients above.".format(row+1, self.leading[row]+1))
                return False

        # no conditions failed
        return True

    def to_reduceret_triangular(self):
        if not self.is_triangular():
            if self.debug:
                print(self)
                print("Has to be in triangular before it can be reduced.")
                
            _ = self.print_log
            self.print_log = False
            self.to_triangular()
            self.print_log = _
            
        if self.debug: print("\n======== Beginning algorithm to get reduced triangular ========")

        # condition 2
        for row in range(self.rows-1, 0-1, -1):
            if self.leading[row] is None: continue
            if not self.lhs[row][self.leading[row]] == 1:
                if self.debug: 
                    print(self)
                    print("The leading coefficient of Row {} is not 1".format(row+1))
                self.leading_to_one(row, self.leading[row])

        # condition 3
        for i, col in enumerate(self.leading[::-1]):
            if col is None: continue
            i = self.rows - i - 1
            for row in range(i-1, -1, -1):
                if self.lhs[row][col] != 0:
                    if self.debug:
                        print(self)
                        print("The leading coefficient at [{}, {}] is not the only one in its column! (Row {}).".format(i+1, col+1, row+1))
                    self.multiply_rows(i, row, col)

        self.log += "\nNow in reduced triangular.\n"
        self.log += self.__str__() + "\n"
        if self.print_log:
            print("\n============================ LOG =================================")
            print(self.log)
        return self


    #### helper functions
    def calculate_leading(self):
        # find the leading cells
        self.leading = []
        for row in range(self.rows):
            self.leading.append(first_nonzero(self.lhs[row]))

    #### row-equivalent operations B E G
    def swap_rows(self, row1, row2):
        # Bytte
        ordering = list(range(self.rows))
        ordering[row1] = row2
        ordering[row2] = row1
        self.lhs = [self.lhs[i] for i in ordering]
        self.rhs = [self.rhs[i] for i in ordering]
        self.log += "\nRow {} <-> Row {}".format(row1+1, row2+1)
        if self.full_log: self.log += self.__str__()
        if self.debug:
            print("The pivotposition is zero, swapping rows.")
            print(self)
        return True

    def multiply_rows(self, row1, row2, pivot_column):
        # Erstatte
        multiplier = self.lhs[row2][pivot_column] / self.lhs[row1][pivot_column]
        sign = "-" if multiplier >= 0 else "+"
        self.lhs[row2] = list(self.lhs[row2][col] - multiplier * self.lhs[row1][col] \
                                for col in range(self.columns))
        self.rhs[row2] = self.rhs[row2] - multiplier * self.rhs[row1]
        
        self.log += "\nRow {} -> Row {} {} {:.3f} * Row {}".format(row2+1,
                                                               row2+1,
                                                               sign,
                                                               multiplier,
                                                               row1+1)
        if self.full_log: self.log += self.__str__() + "\n"
        if self.debug:
            print("Improving row {} by subtracting {} times row {}".format(
                row2+1,
                multiplier,
                row1+1))
            print(self)
        return True

    def leading_to_one(self, row, pivot_column):
        # Gange
        multiplier = 1 / self.lhs[row][pivot_column]
        self.log += "\nRow {} -> Row {} * {:.3f}".format(row+1, row+1, multiplier)
        self.lhs[row] = list(self.lhs[row][col] * multiplier for col in range(self.columns))
        self.rhs[row] = self.rhs[row] * multiplier
        if self.full_log: self.log += self.__str__() + "\n"



# eq1 = Trappe([
#     # lhs
#     [0, 0, 1],
#     [2, 1, 1],
#     [3, 3, 4]],
#     # rhs
#     [0, 1, 10]) \
#     .to_reduceret_trappeform()




# eq131 = Trappe([[2, 3],
#                 [10,9]], [1,11]).to_reduceret_trappeform()

# eq133 = Trappe([[2, -4],
#                 [-1, 5]], [6,0]).to_reduceret_trappeform()


# eq1310 = Trappe([[2, 3, 1],
#                  [4, 7, 5],
#                  [0, -2, 2]], [8, 20, 0], circle_pivots = True).to_reduceret_trappeform()

# eq1311 = Trappe([[2, -3, 0],
#                  [4, -5, 1],
#                  [2, -1, -3]], [3,7,5]).to_reduceret_trappeform()

eq1315 = equation_system([[2, -1, 1],
                          [2, -1, 1],
                          [4, 1, 1]], [0, 0, 2]).to_reduceret_triangular()
