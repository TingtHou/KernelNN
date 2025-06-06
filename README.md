# 📦 KernelNN

KernelNN is an R package for high-dimensional statistical modeling using kernel-based neural networks and variance component estimation. It integrates kernel construction, nonlinear activation, and MINQUE/IMINQUE estimation for both inference and testing.

## 🚀 Key Features
* Build flexible kernel representations of input data
* Apply neural-network-inspired activation to kernels
* Estimate variance components via closed-form MINQUE or iterative MINQUE (IMINQUE)
* Perform hypothesis testing on overall and individual effects
  
## 🔧 Package Usage
### 1. Kernel Construction
Generates a list of kernel matrices from the input data:

```r
KernelList <- KernelGenerator(DataMatrix = Data_Matrix, kernelName = "...")
```
`Data_Matrix`: Input data matrix (features or design).
`kernelName`: The type of kernel matrix to be generated based on the input data matrix. Options include:
* `CAR`: Conditional auto regressive kernel
* `identity`: Identity kernel
* `product`: Product kernel
* `polynomial`: Polynomial kernel
* `Gaussian`: Gaussian kernel
  
### 2.  Apply Activation Function
Transforms the kernel list using a specified activation function:
If the active function is linear:

```r
HiddenKernel <- Linear(KernelList)
```
or polynomial:
```r
HiddenKernel <- polynomial(KernelList)
```
### 3.  Estimate Variance Components
Estimates the variance components using the specified method:
If we use MINQUE0 method
```r
vcs.result <- MINQUE0(KList = HiddenKernel, y = y)
```
or if we use iterative MINQUE method
```r
vcs.result <- IMINQUE(KList = HiddenKernel, y = y, epoch = iteration, threshold = threshold)
```
`epoch`: Maximum number of iterations.
`threshold`: Convergence threshold, 1e-6 commonly.

### 4. Hypothesis Testing

Conduct hypothesis testing on the overall genetic effect as well as selected individual effects. Depending on the chosen variance component estimation method, different testing functions are provided. For more details, please refer to our references.

#### a. If using `MINQUE0`:

```r
Pvalue <- MNQTest0_Chi(y, KList = KList, vcs = vcs.result$vcs, TestingID)
```
`TestingID`: Index or indices of components to test. If `TestingID` is NULL, only the overall effect will be tested.
#### b. If using `IMINQUE`:
```r
Pvalue<-IMNQTest_Normal(y, KList = KList, vcs = vcs.result$vcs, TestingID)
```


### Reference 
Hou, Tingting, Chang Jiang, and Qing Lu. "An Association Test Based on Kernel-Based Neural Networks for Complex Genetic Association Analysis." arXiv preprint arXiv:2312.06669 (2023). 

Hou, Tingting, Chang Jiang, and Qing Lu. "A kernel-based neural network test for high-dimensional sequencing data analysis." arXiv preprint arXiv:2312.02850 (2023).





