goc-solution2-solver Change Log
===============================

### Staged
- Made worker write solution2 files incrementally
- Added automatic deletion of partial solution2 files
- Added correct_voltage_angles! to given solution1
- Improved support for generators with no reactive capability
- Added retries to pmap call in case of node failure
- Cleaned up terminal outputs and worker startup information
- Fixed bug in apply_pg_response!
- Restored solution1 correctness checks

### v0.2.0
- Added multi-node solving support

### v0.1.0
- Initial release

