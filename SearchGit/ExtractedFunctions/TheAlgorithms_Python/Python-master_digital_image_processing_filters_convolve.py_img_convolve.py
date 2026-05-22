def img_convolve(image, filter_kernel):
    height, width = image.shape[0], image.shape[1]
    k_size = filter_kernel.shape[0]
    pad_size = k_size // 2
    # Pads image with the edge values of array.
    image_tmp = pad(image, pad_size, mode="edge")

    # im2col, turn the k_size*k_size pixels into a row and np.vstack all rows
    image_array = im2col(image_tmp, (k_size, k_size))

    #  turn the kernel into shape(k*k, 1)
    kernel_array = ravel(filter_kernel)
    # reshape and get the dst image
    dst = dot(image_array, kernel_array).reshape(height, width)
    return dst