process DECOMPRESS {
    executor 'local'

    input:
    path(compressed)

    output:
    path("${compressed.baseName}")

    script:
    """
    gunzip -c $compressed > ${compressed.baseName}
    """
}
