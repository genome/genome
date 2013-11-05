-- Revert subject.project_part.index_label

BEGIN;

DROP INDEX subject.idx_s_pp_l;

COMMIT;
