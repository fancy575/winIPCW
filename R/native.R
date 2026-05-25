.native_version <- function() {
  .Call("winIPCW_cpp_version", PACKAGE = "winIPCW")
}

.native_pair_ipcw_sums <- function(y_trt, d_trt, y_ctl, d_ctl, g_trt, g_ctl,
                                   tau, eps) {
  .Call(
    "winIPCW_pair_ipcw_sums",
    as.numeric(y_trt),
    as.integer(d_trt),
    as.numeric(y_ctl),
    as.integer(d_ctl),
    as.numeric(g_trt),
    as.numeric(g_ctl),
    as.numeric(tau),
    as.numeric(eps),
    PACKAGE = "winIPCW"
  )
}
