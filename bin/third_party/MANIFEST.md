# Third Party Software Manifest

The distributed binaries are built for the following architectures except where noted:

- Linux/x86_64, RHEL8 or later
- Linux/aarch64, RHEL8 or later
- Apple/x86_64, macOS 10.14 or later
- Apple/arm64, macOS 11.0 or later

Apple binaries (*_Darwin*) are provided in fat/universal form. Artifacts distributed with DAIS-ribosome will use suffixes to denote the appropriate architecture. Copyright and license information is generally copied in the [packaged-citations-licenses](packaged-citations-licenses) folder.

## In DAIS-ribosome and SSWSORT (linked for deployment)

- [GNU Parallel] v20200422
  - Artifacts: `parallel`
  - Requires: system Perl
  - License: [GPL v3]
- [SSW] v1.2.5M
  - Artifacts: `ssw`
  - Custom modifications:
    <https://github.com/sammysheep/Complete-Striped-Smith-Waterman-Library/tree/IRMA%40v1.3>
  - License: [MIT]

[GNU Parallel]: https://www.gnu.org/software/parallel/
[GPL v3]: https://www.gnu.org/licenses/gpl-3.0.txt
[MIT]: https://opensource.org/license/mit
[SSW]: https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
