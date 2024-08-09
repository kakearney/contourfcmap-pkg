[![View Constructor for function_handles on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/45941-constructor-for-function_handles)

[![Donate to Rody](https://i.stack.imgur.com/bneea.png)](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url)

# FEX-function_handle

If you ever found yourself in a situation where you could not or didn't want to add a directory to the MATLAB search path, but still needed speedy access to a function in that directory, this file is for you.
In such cases, having to change directories to gain access to the function is not the best solution: you'd always have to take care to change the current path back to what it was (even on error). Moreover, performance in those cases can be poor if you have to call this function very often.

FUNCTION_HANDLE allows you to create function handles which can successfully be evaluated without loss of performance, even if the function the handle points to is not on the MATLAB search path.

While there often are better ways to accomplish this sort of task (package directories, symbolic links, etc.), there are a few niche cases where these solutions are simply more involved than using this FUNCTION_HANDLE constructor.

Note that FUNCTION_HANDLE overloads a function present in standard MATLAB. This 'native' function is nothing more than documentation (linked to in the help for this constructor) and an error message which says that you cannot use the function to construct handles. As this is exactly what FUNCTION_HANDLE implements, this shadowing is desirable.

Example session:

>> F = function_handle('./path/to/function/myFcn.m')
F =
@myFcn

>> A = function_handle(...
{@cos, '../dir/not/on/path/myFunction.m'})
A =
@cos @myFunction

>> A{1}(pi)
ans =
-1

>> functions(A{1})
ans =
function: 'min'
type: 'simple'
file: ''

>> functions(A{2})
ans =
function: 'myFunction'
type: 'simple'
file: '/fullpath/dir/not/on/path/myFunction.m'

If you like this work, please consider [a donation](https://www.paypal.com/cgi-bin/webscr?cmd=_s-xclick&hosted_button_id=4M7RMVNMKAXXQ&source=url).
