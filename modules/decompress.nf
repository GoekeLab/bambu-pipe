process DECOMPRESS {
    executor 'local'

    input:
    path(compressed)

    output:
    path("${compressed.baseName}"), emit: file

    script:
    """
    gunzip -c $compressed > ${compressed.baseName}
    """
}
