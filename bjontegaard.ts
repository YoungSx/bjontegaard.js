function getBaseLog(x: number, y: number) {
    return Math.log(y) / Math.log(x)
}


function pchip(H: Array<number>, delta: Array<number>) {
    function pchipBoundary(h1: number, h2: number, del1: number, del2: number) {
        let d = ((2 * h1 + h2) * del1 - h1 * del2) / (h1 + h2)

        if (d * del1 < 0)
            d = 0
        else if ((del1 * del2 < 0) && (Math.abs(d) > Math.abs(3 * del1)))
            d = 3 * del1

        return d
    }

    const d: Array<number> = []

    d[0] = pchipBoundary(H[0], H[1], delta[0], delta[1])
    for (let i = 1; i <= H.length; i++) {
        d[i] = (3 * H[i - i] + 3 * H[i]) / ((2 * H[i] + H[i - 1]) / delta[i - 1] + (H[i] + 2 * H[i - 1]) / delta[i])
    }
    const n = H.length
    d[n] = pchipBoundary(H[n - 1], H[n - 2], delta[n - 1], delta[n - 2])

    return d
}


function BDBR(anchor_rate: Array<number>, anchor_distortion: Array<number>, test_rate: Array<number>, test_distortion: Array<number>) {
    function bdRint(rate: Array<number>, distortion: Array<number>, low: number, high: number) {
        const dataLength = rate.length

        // logarithmic
        const logRate: Array<number> = []
        for (let i = 0; i < dataLength; i++)
            logRate[i] = getBaseLog(10, rate[i])

        const logDistortion: Array<number> = []
        for (let i = 0; i < dataLength; i++)
            logDistortion[i] = distortion[i]

        const H: Array<number> = []
        const delta: Array<number> = []
        for (let i = 0; i < dataLength - 1; i++) {
            H[i] = logDistortion[i + 1] - logDistortion[i]
            delta[i] = (logRate[i + 1] - logRate[i]) / H[i] // gradient of sample point
        }

        const d: Array<number> = pchip(H, delta)

        const c: Array<number> = []
        const b: Array<number> = []
        for (let i = 0; i < dataLength - 1; i++) {
            c[i] = (2 * delta[i] - 2 * d[i] - d[i + 1]) / H[i]
            b[i] = (d[i] - 2 * delta[i] + d[i + 1]) / (Math.pow(H[i], 2))
        }

        let s0: number // starting point of the line segment
        let s1: number // ending point of the line segment
        let integral = 0

        for (let i = 0; i < dataLength - 1; i++) {
            s0 = logDistortion[i]
            s1 = logDistortion[i + 1]

            // normalized to within range
            s0 = Math.max(s0, low)
            s0 = Math.min(s0, high)
            s1 = Math.max(s1, low)
            s1 = Math.min(s1, high)

            // move to the coordinate axis
            s0 = s0 - logDistortion[i]
            s1 = s1 - logDistortion[i]

            if (s1 > s0) {
                integral = integral + (s1 - s0) * logRate[i] +
                    (Math.pow(s1, 2) - Math.pow(s0, 2)) * d[i] / 2 +
                    (Math.pow(s1, 3) - Math.pow(s0, 3)) * c[i] / 3 +
                    (Math.pow(s1, 4) - Math.pow(s0, 4)) * b[i] / 4
            }
        }

        return integral
    }

    const distortionMin = Math.max(Math.min(...anchor_distortion), Math.min(...test_distortion))
    const distortionMax = Math.min(Math.max(...anchor_distortion), Math.max(...test_distortion))

    const integralAnchor = bdRint(anchor_rate, anchor_distortion, distortionMin, distortionMax)
    const integralTest = bdRint(test_rate, test_distortion, distortionMin, distortionMax)

    const avg = (integralTest - integralAnchor) / (distortionMax - distortionMin)

    return Math.pow(10, avg) - 1
}


function BDPSNR(anchor_rate: Array<number>, anchor_distortion: Array<number>, test_rate: Array<number>, test_distortion: Array<number>) {
    function bdRint(rate: Array<number>, distortion: Array<number>, low: number, high: number) {
        const dataLength = rate.length

        // logarithmic
        const logRate: Array<number> = []
        for (let i = 0; i < dataLength; i++)
            logRate[i] = getBaseLog(10, rate[i])

        const logDistortion: Array<number> = []
        for (let i = 0; i < dataLength; i++)
            logDistortion[i] = distortion[i]

        const H: Array<number> = []
        const delta: Array<number> = []
        for (let i = 0; i < dataLength - 1; i++) {
            H[i] = logRate[i + 1] - logRate[i]
            delta[i] = (logDistortion[i + 1] - logDistortion[i]) / H[i] // gradient between 2 points
        }

        const d: Array<number> = pchip(H, delta) // gradient of sample point

        const c: Array<number> = []
        const b: Array<number> = []
        for (let i = 0; i < dataLength - 1; i++) {
            c[i] = (2 * delta[i] - 2 * d[i] - d[i + 1]) / H[i]
            b[i] = (d[i] - 2 * delta[i] + d[i + 1]) / (Math.pow(H[i], 2))
        }

        let s0: number // starting point of the line segment
        let s1: number // ending point of the line segment
        let integral = 0

        for (let i = 0; i < dataLength - 1; i++) {
            s0 = logRate[i]
            s1 = logRate[i + 1]

            // normalized to within range
            s0 = Math.max(s0, low)
            s0 = Math.min(s0, high)
            s1 = Math.max(s1, low)
            s1 = Math.min(s1, high)

            // move to the coordinate axis
            s0 = s0 - logRate[i]
            s1 = s1 - logRate[i]

            if (s1 > s0) {
                integral = integral + (s1 - s0) * logDistortion[i] +
                    (Math.pow(s1, 2) - Math.pow(s0, 2)) * d[i] / 2 +
                    (Math.pow(s1, 3) - Math.pow(s0, 3)) * c[i] / 3 +
                    (Math.pow(s1, 4) - Math.pow(s0, 4)) * b[i] / 4
            }
        }

        return integral
    }

    const rateMin = Math.max(Math.min(...anchor_rate), Math.min(...test_rate))
    const rateMax = Math.min(Math.max(...anchor_rate), Math.max(...test_rate))

    const logRateMin = getBaseLog(10, rateMin)
    const logRateMax = getBaseLog(10, rateMax)

    const integralAnchor = bdRint(anchor_rate, anchor_distortion, logRateMin, logRateMax)
    const integralTest = bdRint(test_rate, test_distortion, logRateMin, logRateMax)

    const avg = (integralTest - integralAnchor) / (logRateMax - logRateMin)

    return avg
}


function test() {
    const rate1 = [1175629, 1872841, 3068206, 4824244]
    const rate2 = [1159665, 1884065, 3023406, 4684530]

    const psnr1 = [39.077549, 40.778417, 42.63382, 44.351552]
    const psnr2 = [39.245824, 40.996712, 42.753412, 44.414841]

    console.log(`BDBR: ${BDBR(rate1, psnr1, rate2, psnr2)}`)
    console.log(`BDPSNR: ${BDPSNR(rate1, psnr1, rate2, psnr2)}`)
}

test()
